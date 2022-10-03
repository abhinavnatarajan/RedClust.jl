using StatsBase: autocor, wsample, levelsmap
using Distributions: Dirichlet, MvNormal, pdf
using LinearAlgebra: I
using Random: rand

"""
    adjacencymatrix(clusts::ClustLabelVector)
Returns the `n`x`n` adjacency matrix corresponding to the given cluster label vector `clusts`, where `n = length(clusts)`. 
""" 
function adjacencymatrix(
    clusts::ClustLabelVector
    )::Matrix{Bool}
	clusts .== clusts'
end



function uppertriangle(
    M::Matrix{T}
    )::Vector{T} where {T}
	n = size(M, 1)
	[M[i, j] for i in 1:n for j in (i + 1):n]
end

"""
    sortlabels(x::ClustLabelVector)
Returns a cluster label vector `y` such that `x` and `y` have the same adjacency structure and labels in `y` occur in sorted ascending order.
"""
function sortlabels(
    x::ClustLabelVector
    )::ClustLabelVector 
    temp = levelsmap(x)
    [temp[i] for i in x]
end

function numpairs(n::Int)::Int
    Integer(n * (n - 1) / 2)
end

function sample_logweights(logprobs::Vector{Float64})::Int
	logprobs = logprobs .- minimum(logprobs)
    u = rand(Uniform(), length(logprobs))
    argmax(-log.(-log.(u)) + logprobs)
end

# Compute the integrated autocorrelation coefficient
function iac_ess_acf(x::AbstractVector{<:Real})
    acf = autocor(x)
    iac = sum(acf) * 2
    ess = length(x) / iac
    return (iac = iac, ess = ess, acf = acf)
end

"""
    generatemixture(N, K; α = K, dim = K, radius = 1, σ = 0.1)

Generates a multivariate Normal mixture, with kernel weights generated from a Dirichlet prior. The kernels are centred at the vertices of a `dim`-dimensional simplex with edge length `radius`.

# Arguments
- `N::Integer`: number of observations to generate.
- `K::Integer`: number of mixture kernels.
- `α::Float64 = K`: parameter for the Dirichlet prior.
- `dim::Integer = K`: dimension of the observations.
- `radius::Float64 = 1`: radius of the simplex whose vertices are the kernel means.
- `σ::Float64 = 0.1 `: variance of each kernel.

# Returns

Named tuple containing the following fields-

- `pnts::Vector{Vector{Float64}}`: a vector of `N` observations.
- `distM::Matrix{Float64}`: an `N`x`N` matrix of pairwise Euclidean distances between the observations.
- `clusts::ClustLabelVector`: vector of `N` cluster assignments.
- `probs::Float64`: vector of `K` cluster weights generated from the Dirichlet prior, used to generate the observations.
- `oracle_coclustering::Matrix{Float64}`: `N`x`N` matrix of co-clustering probabilities, calculated assuming full knowledge of the cluster centres and cluster weights.
"""
function generatemixture(N, K; α = K, dim = K, radius = 1, σ = 0.1) 
    probs = rand(Dirichlet(K, α))
    clusts = sort(wsample(1:K, probs, N))

	# generate cluster centres
	clust_centres = fill(zeros(dim), K)
	for i = 1:K
		clust_centres[i] = cat(zeros(i - 1), radius, zeros(dim - i); dims = 1)
    end

	# generate points
    Σ = σ^2 * I(dim)
	pnts = fill(zeros(dim), N)
	oracle_posterior = zeros(K, N)
	for i = 1:N
		pnts[i] = rand(MvNormal(clust_centres[clusts[i]], Σ))
	end
    # calculate oracle
    numiters = 1000
    oracle_posterior = zeros(K, N)
    oracle_coclustering = zeros(N, N)
    for c in 1:numiters
        tempprobs = rand(Dirichlet(K, α))
        for i = 1:N
            for j = 1:K
                oracle_posterior[j, i] = tempprobs[j] * pdf(MvNormal(clust_centres[j], Σ), pnts[i])
            end
            oracle_posterior[:, i] .= oracle_posterior[:, i] ./ sum(oracle_posterior[:, i])
        end
        oracle_coclustering += oracle_posterior' * oracle_posterior
    end
    oracle_coclustering ./= numiters #oracle coclustering matrix
	distM = [pnts[i][j] for i in 1:N, j in 1:dim]' |> 
    (x -> pairwise(Euclidean(), x, dims = 2))
    return (pnts = pnts, distM = distM, clusts = clusts, probs = probs, oracle_coclustering = oracle_coclustering)
end

function makematrix(x::Vector{Vector{T}})::Matrix where {T}
    [x[i][j] for j in 1:length(x[1]), i in 1:length(x)]
end



