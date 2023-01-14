using Clustering: kmeans, kmedoids
using Distances: pairwise, Euclidean
using Distributions: fit_mle, Gamma, shape, rate, Distributions, Beta
using ProgressBars: ProgressBar
using StatsBase: counts, std
using SpecialFunctions: logbeta

"""
    fitprior(data, algo, diss = false; 
	Kmin = 1, 
	Kmax = Int(floor(size(data)[end] / 2), 
	verbose = true)

Determines the best prior hyperparameters from the data. A notional clustering is obtained using k-means or k-medoids and the elbow method, and the distances are split into within-cluster distances and inter-cluster distances based on the notional clustering. These distances are then used to fit the prior hyperparameters using MLE and empirical Bayes sampling.   

# Required Arguments
- `data::Union{Vector{Vector{Float64}}, Matrix{Float64}}`: can either be a vector of (possibly multi-dimensional) observations, or a matrix with each column an observation, or a square matrix of pairwise dissimilarities. 
- `algo::String`: must be one of `"k-means"` or `"k-medoids"`.

# Optional Arguments
- `diss::bool`: if true, `data` will be assumed to be a pairwise dissimilarity matrix. 
- `Kmin::Integer`: minimum number of clusters to scan for the elbow method.  
- `Kmax::Integer`: maximum number of clusters to scan for the elbow method. If left unspecified, it is set to half the number of observations.
- `verbose::Bool`: if false, disables all info messages and progress bars. 

# Returns
An object of type [`PriorHyperparamsList`](@ref).
"""
function fitprior(
	data::Union{Vector{Vector{Float64}}, Matrix{Float64}},
	algo::String, 
	diss::Bool = false; 
	Kmin::Integer = 1, 
	Kmax::Integer = Int(floor(size(data)[end] / 2)),
	useR::Bool = false, 
	verbose::Bool = true 
)
	ostream = verbose ? stdout : devnull
	printstyled(ostream, "Fitting prior hyperparameters\n"; bold = true)
	if data isa Vector{Vector{Float64}}
		if diss
			throw(ArgumentError("diss = true but data is not a dissimilarity matrix. Assuming that the data is a vector of observations."))
		end
		x = makematrix(data)
		diss = false
	else
		x = data
	end
	N = size(x, 2)

	diss ? println(ostream, "Input: pairwise dissimilarities between $N observations.") :
		println(ostream, "Input: $N observations of dimension $(size(x, 1)).")

	# Input validation
	diss && size(x, 1) != size(x, 2) && throw(ArgumentError("Supplied dissimilarity matrix is not square."))
	algo == "k-means" && diss && throw(ArgumentError("Cannot use algorithm `k-means` with a dissimilarity matrix."))
	algo != "k-means" && algo != "k-medoids" && throw(ArgumentError("Algo must be 'k-means' or 'k-medoids'."))
	!(1 ≤ Kmin && Kmin ≤ Kmax && Kmax ≤ N) && throw(ArgumentError("Kmin and Kmax must satisfy 1 ≤ Kmin ≤ Kmax ≤ N"))
	dissM = diss ? x : pairwise(Euclidean(), x, dims = 2)

	# Get notional clustering
	println(ostream, "Computing notional clustering.")
	objective = Vector{Float64}(undef, Kmax - Kmin + 1)
	if algo == "k-means"
		clustfn = kmeans
		input = x
	elseif algo == "k-medoids"
		input = dissM
		clustfn = kmedoids
	end

	@inbounds for k in ProgressBar(1:(Kmax-Kmin+1), output_stream = ostream)
		temp = clustfn(input, k; maxiter=1000)
		objective[k] = temp.totalcost
		!temp.converged && @warn "Clustering did not converge at K = $k"
	end
	elbow = detectknee(Kmin:Kmax, objective)[1]
	K = elbow
	notionalclustering = clustfn(input, K; maxiter=1000).assignments
	notionaladjmatrix = adjacencymatrix(notionalclustering)
	A = uppertriangle(dissM)[uppertriangle(notionaladjmatrix) .== 1]
	B = uppertriangle(dissM)[uppertriangle(notionaladjmatrix) .== 0]

	# Compute partition prior parameters
	println(ostream, "Computing partition prior hyperparameters.")
	clustsizes = counts(notionalclustering)
	temp = sample_rp(clustsizes; verbose = verbose)
	proposalsd_r = std(temp.r)
	fitr = fit_mle(Gamma, temp.r)
	η = shape(fitr)
	σ = rate(fitr)
	fitp = fit_mle(Beta, temp.p)
	u, v = Distributions.params(fitp)

	# Compute likelihood parameters
	println(ostream, "Computing likelihood hyperparameters.")
	if K == N # A is empty
		@warn "Got a notional clustering of entirely singletons. Falling back to defaults for cohesion parameters."
		δ1 = 1
		α = 1
		β = 1
	else
		fitA = fit_mle(Gamma, A)
		δ1 = shape(fitA)
		α = length(A) * δ1
		β = sum(A)
	end
	if K == 1
		@warn "Got a notional clustering with a single cluster. Falling back to defaults for repulsion parameters."
		δ2 = 1
		ζ = 1
		γ = 1
	else
		fitB = fit_mle(Gamma, B)
		δ2 = shape(fitB) 
		ζ = length(B) * δ2
		γ = sum(B)
	end

	params = PriorHyperparamsList(
		δ1 = δ1,
		δ2 = δ2,
		α = α,
		β = β,
		ζ = ζ,
		γ = γ, 
		η = η,
		σ = σ,
		proposalsd_r = proposalsd_r,
		u = u, 
		v = v,
		K_initial = K
	)
	return params
end

"""
	sampledist(params::PriorHyperparamsList, type::String, numsamples = 1)::Vector{Float64}

Generate a vector of samples of length `numsamples` from the prior predictive distribution on the distances, as encapsulated in `params`. `type` must be either `"intercluster"` or `"intracluster"`.
"""
function sampledist(params::PriorHyperparamsList, type::String, numsamples::Int = 1)::Vector{Float64}
	# Input validation
	if type ∉ ["intercluster", "intracluster"]
		throw(ArgumentError("type must be either \"intercluster\" or \"intracluster\"."))
	end
	if numsamples < 1
		throw(ArgumentError("numsamples must be a positive integer."))
	end
	if type == "intracluster"
		a = params.α
		b = params.β
		δ = params.δ1
	else
		a = params.ζ
		b = params.γ
		δ = params.δ2
	end
	dist1 = Gamma(a, 1/b)
	samples = Vector{Float64}(undef, numsamples)
	for i = 1:numsamples
		dist2 = Gamma(δ, 1/rand(dist1))
		samples[i] = rand(dist2)
	end
	samples
end

"""
	sampleK(params::PriorHyperparamsList, numsamples::Int, n::Int)::Vector{Int}

Returns a vector of length `numsamples` containing samples of ``K`` (number of clusters) generated from its marginal prior predictive distribution inferred from `params`. The parameter `n` is the number of observations in the model.
"""
function sampleK(params::PriorHyperparamsList, numsamples::Int, n::Int)::Vector{Int}
	# Input validation
	if n < 1 
		throw(ArgumentError("n must be a positive integer."))
	end
	if numsamples < 1
		throw(ArgumentError("numsamples must be a positive integer."))
	end
	samples = Vector{Int}(undef, numsamples)
	rdist = Gamma(params.η, 1/params.σ)
	pdist = Beta(params.u, params.v)
	K = 1:(n-1)
	logprobs = zeros(n)
	@inbounds for i = 1:numsamples
		r = rand(rdist)
		p = rand(pdist)
		@. logprobs[1:(n-1)] = (r*K) * log(1-p) + (n-K)* log(p) - log(n-K) - logbeta(r*K, n-K)
		logprobs[n] = r*n * log(1-p)
		samples[i] = sample_logweights(logprobs)
	end
	return samples
end

function detectknee(
	xvalues::AbstractVector{<:Real}, 
	yvalues::AbstractVector{<:Real}
	)::NTuple{2,Real}
	# Max values to create line
	ind = sortperm(xvalues)
	x = xvalues[ind]
	y = yvalues[ind]
	x1 = x[1]; x2 = x[end]
	y1 = y[1]; y2 = y[end]

	# Line between extreme values
	a = (y2-y1)/(x2-x1)
	b = y1 - a*x1

	# Distance from points to line
	distances = [abs(a * x[i] + b - y[i]) / sqrt(a^2 + 1) for i in eachindex(x)]
	# Max distance point
	maxind = argmax(distances)
	x[maxind], y[maxind]
end