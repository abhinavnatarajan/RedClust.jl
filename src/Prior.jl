using Distributions: fit_mle, Gamma, shape, rate, Distributions
using Clustering: kmeans, kmedoids
using Distances: pairwise, Euclidean
using StatsBase: counts, std
using ProgressBars: ProgressBar
using RCall: rcopy, @R_str

"""
    fitprior(data, algo, dist = false, Kmin = 1, Kmax = Int(floor(size(data)[end] / 2))

Compute the prior hyperparameters from the data. 

# Arguments
- `data::Union{Vector{Vector{Float64}}, Matrix{Float64}}`: can either be a vector of (possibly multi-dimensional) observations, or a matrix with each column an observation, or a square matrix of pairwise dissimilarities. 
- `algo::String`: must be one of `"k-means"` or `"k-medoids"`.
- `dist::bool = false`: if true, `data` will be assumed to be a pairwise dissimilarity matrix. 
- `Kmin::Integer`: minimum number of clusters.
- `Kmax::Integer = Int(floor(size(data)[end] / 2))`: maximum number of clusters. If left unspecified, it is set to half the number of observations.

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
	dist::Bool = false, 
	Kmin = 1, 
	Kmax = Int(floor(size(data)[end] / 2)) 
)
	if data isa Vector{Vector{Float64}}
		x = makematrix(data)
		dist = false
	else
		x = data
	end
	N = size(x, 2)

	# Input validation
	if dist && size(x, 1) != size(x, 2)
		error("Supplied distance matrix is not square.")
	end
	if algo == "k-means" && dist
		error("Cannot use algorithm `k-means` with a distance matrix.")
	end

	if !dist
		distM = pairwise(Euclidean(), x, dims = 2)
	elseif dist
		distM = x
	end
	

	# Get notional clustering
	print("Computing notional clustering.\n")
	objective = Vector{Float64}(undef, Kmax - Kmin + 1)
	if algo == "k-means"
		if dist
			error("Cannot use algo k-means with distance matrix.")
		end
		temp = x'
		objective = rcopy(
			R"""
			sapply($Kmin:$Kmax, 
			function(k){kmeans($temp, centers = k, nstart = 20)$tot.withinss})
			"""
		) 
	elseif algo == "k-medoids"
		objective = rcopy(
			R"""
			library(cluster)
			sapply($Kmin:$Kmax, 
			function(k){pam($distM, k = k, diss = T, pamonce = 5)$objective[["swap"]]})
			"""
		)
	else
		error("Algo must be 'k-means' or 'k-medoids'.")
	end

	# for k in ProgressBar(1:(Kmax-Kmin+1))
	# 	temp = clustfn(input, k; maxiter=1000)
	# 	objective[k] = temp.totalcost
	# 	if !temp.converged
	# 		@warn "Clustering did not converge at K = $k"
	# 	end
	# end
	elbow = detectknee(Kmin:Kmax, objective)[1]
	K = elbow
	if (algo == "k-means")
		notionalclustering = rcopy(
			R"""
			kmeans($temp, centers = $K, nstart = 20)
			"""
		)[:cluster]
	elseif algo == "k-medoids"
		notionalclustering = rcopy(
			R"""
			library(cluster)
			pam($distM, k = $K, diss = T, pamonce = 5)
			"""
		)[:clustering]
	end
	# notionalclustering = clustfn(input, K; maxiter=1000)
	# notionalclustering = notionalclustering.assignments
	notionaladjmatrix = adjacencymatrix(notionalclustering)
	A = uppertriangle(distM)[uppertriangle(notionaladjmatrix) .== 1]
	B = uppertriangle(distM)[uppertriangle(notionaladjmatrix) .== 0]

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