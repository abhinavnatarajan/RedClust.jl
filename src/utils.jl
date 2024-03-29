# use the Gumbel-max trick to sample from a vector of discrete log-probabilities
@inline function sample_logweights(logprobs::AbstractVector{Float64})::Int
    logprobs .-= minimum(logprobs)
    u = rand(Uniform(), length(logprobs))
    argmax(-log.(-log.(u)) .+ logprobs)
end

# Fast versions of vector and matrix sums
@inline function matsum(x::AbstractMatrix{Float64}, inds1::AbstractVector{Int}, inds2::AbstractVector{Int})
    ans = 0.0
    @turbo for i in eachindex(inds1)
        for j in eachindex(inds2)
            ans += x[inds1[i], inds2[j]]
        end
    end
    ans
end
@inline function matsum(x::AbstractMatrix{Float64})
    ans = 0.0
    @turbo for i in eachindex(x)
        ans += x[i]
    end
    ans
end
@inline function vecsum(x::AbstractVector{Float64}, inds::AbstractVector{Int})
    ans = 0.0
    @turbo for i in eachindex(inds)
        ans += x[inds[i]]
    end
    ans
end
@inline function vecsum(x::AbstractVector{Float64})
    ans = 0.0
    @turbo for i in eachindex(x)
        ans += x[i]
    end
    ans
end

# Compute the integrated autocorrelation coefficient
function iac_ess_acf(x::AbstractVector{<:Real})
    acf = autocor(x)
    iac = sum(acf) * 2
    ess = length(x) / iac
    return (iac=iac, ess=ess, acf=acf)
end

function uppertriangle(
    M::Matrix{T}
)::Vector{T} where {T}
    n = size(M, 1)
    [M[i, j] for i in 1:n for j in (i+1):n]
end

"""
adjacencymatrix(clusts::ClustLabelVector) -> Matrix{Bool}
Returns the `n`×`n` adjacency matrix corresponding to the given cluster label vector `clusts`, where `n = length(clusts)`.
"""
function adjacencymatrix(
    clusts::ClustLabelVector
)::Matrix{Bool}
    clusts .== clusts'
end

"""
sortlabels(x::ClustLabelVector) -> ClustLabelVector
Returns a cluster label vector `y` such that `x` and `y` have the same adjacency structure and labels in `y` occur in sorted ascending order.
"""
function sortlabels(
    x::ClustLabelVector
)::ClustLabelVector
    temp = levelsmap(x)
    [temp[i] for i in x]
end

"""
generatemixture(N, K; [α, dim, radius, σ, rng])

Generates a multivariate Normal mixture, with kernel weights generated from a Dirichlet prior. The kernels are centred at the vertices of a `dim`-dimensional simplex with edge length `radius`.

# Required Arguments
- `N::Integer`: number of observations to generate.
- `K::Integer`: number of mixture kernels.

# Optional Arguments
- `α::Float64 = K`: parameter for the Dirichlet prior.
- `dim::Integer = K`: dimension of the observations.
- `radius::Float64 = 1`: radius of the simplex whose vertices are the kernel means.
- `σ::Float64 = 0.1 `: variance of each kernel.
- `rng`: a random number generator to use, or an integer to seed the default random number generator with. If not provided, the default RNG provided by the Random.jl package will be used with default seeding.

# Returns
Named tuple containing the following fields-

- `points::Vector{Vector{Float64}}`: a vector of `N` observations.
- `distancematrix::Matrix{Float64}`: an `N`×`N` matrix of pairwise Euclidean distances between the observations.
- `clusts::ClustLabelVector`: vector of `N` cluster assignments.
- `probs::Float64`: vector of `K` cluster weights generated from the Dirichlet prior, used to generate the observations.
- `oracle_coclustering::Matrix{Float64}`: `N`×`N` matrix of co-clustering probabilities, calculated assuming full knowledge of the cluster centres and cluster weights.
"""
function generatemixture(N::Integer, K::Integer; α::Real=K, dim::Real=K, radius::Real=1, σ::Real=0.1, rng::Union{AbstractRNG,<:Integer}=TaskLocalRNG())
    # input validation
    N < 1 && throw(ArgumentError("N must be greater than 1."))
    (K < 1 || K > N) && throw(ArgumentError("K must satisfy 1 ≤ K ≤ N."))
    α ≤ 0 && throw(ArgumentError("α must be positive."))
    dim < K && throw(ArgumentError("dim must be ≥ K."))
    radius ≤ 0 && throw(ArgumentError("radius must be positive."))
    σ ≤ 0 && throw(ArgumentError("σ must be positive."))
    if rng isa Int
        rng = seed!(rng)
    end

    probs = rand(rng, Dirichlet(K, α))
    clusts = sort(wsample(rng, 1:K, probs, N))

    # generate cluster centres
    clust_centres = fill(zeros(dim), K)
    for i = 1:K
        clust_centres[i] = cat(zeros(i - 1), radius, zeros(dim - i); dims=1)
    end

    # generate points
    Σ = σ^2 * I(dim)
    pnts = fill(zeros(dim), N)
    oracle_posterior = zeros(K, N)
    for i = 1:N
        pnts[i] = rand(rng, MvNormal(clust_centres[clusts[i]], Σ))
    end
    # calculate oracle
    numiters = 5000
    oracle_posterior = zeros(K, N)
    oracle_coclustering = zeros(N, N)
    for _ in 1:numiters
        tempprobs = rand(rng, Dirichlet(K, α))
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
            (x -> pairwise(Euclidean(), x, dims=2))
    return (points=pnts, distancematrix=distM, clusts=clusts, probs=probs, oracle_coclustering=oracle_coclustering)
end

"""
makematrix(x::AbstractVector{<:AbstractVector}) -> Matrix

Convert a vector of vectors into a matrix, where each vector becomes a column in the matrix.
"""
function makematrix(x::AbstractVector{<:AbstractVector})::Matrix
    [x[i][j] for j in 1:length(x[1]), i in 1:length(x)]
end

function prettytime(t)
    if t < 1e-6 # nanoseconds
        out = @sprintf "%.2f ns" t * 1e9
    elseif t ≥ 1e-6 && t < 1e-3 # microseconds
        out = @sprintf "%.2f μs" t * 1e6
    elseif t ≥ 1e-3 && t < 1 # milliseconds
        out = @sprintf "%.2f ms" t * 1e3
    elseif t < 60 # seconds
        out = @sprintf "%.2f s" t
    else
        x = Dates.canonicalize(Dates.Second(Integer(floor(t))))
        out = ""
        for i in 1:length(x.periods)
            out *= string(x.periods[i].value) * " "
            if typeof(x.periods[i]) == Second
                out *= "s"
            elseif typeof(x.periods[i]) == Minute
                out *= "min" * (x.periods[i].value != 1 ? "s" : "")
            elseif typeof(x.periods[i]) == Hour
                out *= "hr" * (x.periods[i].value != 1 ? "s" : "")
            elseif typeof(x.periods[i]) == Day
                out *= "day" * (x.periods[i].value != 1 ? "s" : "")
            end
            if i != length(x.periods)
                out *= " "
            end
        end
    end
    out
end

prettynumber(x) = x < 1 ? @sprintf("%.3e", x) : @sprintf("%.3f", x)
