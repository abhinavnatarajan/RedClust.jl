using Distributions: fit_mle, Gamma, shape, rate, Distributions
using Clustering: kmeans, kmedoids
using Distances: pairwise, Euclidean
using StatsBase: counts, std
using ProgressBars: ProgressBar
using RCall: rcopy, @R_str

"""
    fitprior(data, algo; diss = false, Kmin = 1, Kmax = Int(floor(size(data)[end] / 2), useR = true)

Determines the best prior hyperparameters from the data. A notional clustering is obtained using k-means or k-medoids, and the distances are split into within-cluster distances and inter-cluster distances based on the notional clustering. These distances are then used to fit the prior hyperparameters using MLE and empirical Bayes sampling.   

# Arguments
- `data::Union{Vector{Vector{Float64}}, Matrix{Float64}}`: can either be a vector of (possibly multi-dimensional) observations, or a matrix with each column an observation, or a square matrix of pairwise dissimilarities. 
- `algo::String`: must be one of `"k-means"` or `"k-medoids"`.
- `diss::bool = false`: if true, `data` will be assumed to be a pairwise dissimilarity matrix. 
- `Kmin::Integer`: minimum number of clusters.
- `Kmax::Integer = Int(floor(size(data)[end] / 2))`: maximum number of clusters. If left unspecified, it is set to half the number of observations.
- `useR::Bool = false`: if `false`, will use the `kmeans` or `kmedoids` from the Julia package [`Clustering.jl`](https://juliastats.org/Clustering.jl/stable/). If `true`, will use `kmeans` or `cluster::pam` in R. 

# Returns
A named tuple containing the following
- `params::PriorHyperparamsList`: the fitted hyperparameters.
- `K::Integer`: the number of clusters in the notional clustering.
- `notionalclustering::ClustLabelVector`: the cluster labels in the notional clustering.

# See also
[`PriorHyperparamsList`](@ref)
"""
function fitprior(
	data::Union{Vector{Vector{Float64}}, Matrix{Float64}},
	algo::String, 
	diss::Bool = false; 
	Kmin = 1, 
	Kmax = Int(floor(size(data)[end] / 2)),
	useR = false 
)
	if data isa Vector{Vector{Float64}}
		if diss
			@warn "diss = true but data is not a dissimilarity matrix. Assuming that the data is a vector of observations."
		end
		x = makematrix(data)
		diss = false
	else
		x = data
	end
	N = size(x, 2)

	if !diss 
		print("Input: $N observations.\n")
	else
		print("Input: pairwise dissimilarities of $N observations.\n")
	end

	# Input validation
	if diss && size(x, 1) != size(x, 2)
		throw(ArgumentError("Supplied dissimilarity matrix is not square."))
	end
	if algo == "k-means" && diss
		throw(ArgumentError("Cannot use algorithm `k-means` with a dissimilarity matrix."))
	end
	if algo != "k-means" && algo != "k-medoids"
		throw(ArgumentError("Algo must be 'k-means' or 'k-medoids'."))
	end
	if Kmin < 1
		@warn ("Kmin must be positive, setting Kmin = 1.")
		Kmin = 1
	end
	if Kmax > N
		@warn "Kmax is larger than the number of points (did you accidentally tranpose the input?). Setting Kmax = number of points."
		Kmax = N
	end

	if !diss
		dissM = pairwise(Euclidean(), x, dims = 2)
	else
		dissM = x
	end
	

	# Get notional clustering
	print("Computing notional clustering.\n")
	objective = Vector{Float64}(undef, Kmax - Kmin + 1)
	if algo == "k-means"
		if useR
			input = x'
			objective = rcopy(
				R"""
				sapply($Kmin:$Kmax, 
				function(k){kmeans($input, centers = k, nstart = 20)$tot.withinss})
				"""
				) 
		else
			clustfn = kmeans
			input = x
		end
	elseif algo == "k-medoids"
		input = dissM
		if useR
			objective = rcopy(
				R"""
				library(cluster)
				sapply($Kmin:$Kmax, 
				function(k){pam($dissM, k = k, diss = T, pamonce = 5)$objective[["swap"]]})
				"""
				)
		else
			clustfn = kmedoids
		end
	end

	if !useR
		for k in ProgressBar(1:(Kmax-Kmin+1))
			temp = clustfn(input, k; maxiter=1000)
			objective[k] = temp.totalcost
			# if !temp.converged
			# 	@warn "Clustering did not converge at K = $k"
			# end
		end
	end
	elbow = detectknee(Kmin:Kmax, objective)[1]
	K = elbow
	if useR
		if (algo == "k-means")
			notionalclustering = rcopy(
				R"""
				kmeans($input, centers = $K, nstart = 20)
				"""
			)[:cluster]
		elseif algo == "k-medoids"
			notionalclustering = rcopy(
				R"""
				library(cluster)
				pam($input, k = $K, diss = T, pamonce = 5)
				"""
			)[:clustering]
		end
	else
		notionalclustering = clustfn(input, K; maxiter=1000).assignments
	end
	notionaladjmatrix = adjacencymatrix(notionalclustering)
	A = uppertriangle(dissM)[uppertriangle(notionaladjmatrix) .== 1]
	B = uppertriangle(dissM)[uppertriangle(notionaladjmatrix) .== 0]

	# Compute likelihood parameters
	print("Computing likelihood hyperparameters.\n")
	fitA = fit_mle(Gamma, A)
	fitB = fit_mle(Gamma, B)
	δ1 = shape(fitA)
	α = length(A) * δ1
	β = sum(A)
	δ2 = shape(fitB) 
	ζ = length(B) * δ2
	γ = sum(B)

	# Compute partition prior parameters
	print("Computing partition prior hyperparameters.\n")
	clustsizes = counts(notionalclustering)
	temp = sample_rp(clustsizes)
	proposalsd_r = std(temp.r)
	fitr = fit_mle(Gamma, temp.r)
	η = shape(fitr)
	σ = rate(fitr)
	fitp = fit_mle(Beta, temp.p)
	u, v = Distributions.params(fitp)

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
	@NamedTuple{
		params, K, notionalclustering, 
		objective}(
			(params, K, notionalclustering, objective)
			)
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