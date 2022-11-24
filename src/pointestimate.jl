using Clustering: randindex, mutualinfo, varinfo
using StatsBase: entropy

"""
    getpointestimate(samples::MCMCResult; method = "MAP", loss = "VI")

Computes a point estimate from a vector of samples of cluster allocations by searching for a minimiser of the posterior expectation of some loss function. 

# Arguments 
- `method::String `: must be one of the following -
    - `"MAP"`: maximum a posteriori.
    - `"MLE"`: maximum likelihood estimation.
    - `"MPEL"`: minimum posterior expected loss.
    `"MAP"` and `"MLE"` search among the MCMC samples for the clustering with the maximum log posterior and log likelihood respectively. `"MPEL"` searches for a clustering that minimises the posterior expected loss of some loss function specified by the `loss` argument. The search space is the set of samples in `samples`. 
- `loss`: Determines the loss function used for the `"MPEL"` method. Must be either be a string or a function. If specified as a string, must be one of `"binder"` (Binder loss), `"omARI"` (one minus the Adjusted Rand Index), `"VI"` (Variation of Information distance), or `"ID"` (Information Distance). If specified as a function, must have a method defined for `(x::Vector{Int}, y::Vector{Int}) -> Real`. 

# Returns
Returns a tuple `(clust, i)` where `clust` is a clustering in `samples` and `i` is its sample index. 
"""
function getpointestimate(samples::MCMCResult; method::String = "MAP", loss::Union{String, Function} = "VI")::Tuple{ClustLabelVector, Int}
    # input validation
    if method == "MPEL" && typeof(loss) == String
        if loss ∉ ["binder", "omARI", "VI", "ID"]
            throw(ArgumentError("Invalid loss function specifier."))
        end
    end
    if method ∉ ["MAP", "MLE", "MPEL"]
        throw(ArgumentError("Invalid method specifier."))
    end

    if method == "MAP" # maximum a posteriori
        i = argmax(samples.logposterior)
        return (samples.clusts[i], i)
    elseif method == "MLE" # maximum likelihood estimation
        i = argmax(samples.loglik)
        return (samples.clusts[i], i)
    elseif method == "MPEL" # minimum posterior expectated loss
        # Set loss function
        if loss isa Function
            lossfn = loss
        elseif typeof(loss) == String
            if loss == "binder" 
                lossfn = (x::ClustLabelVector, y::ClustLabelVector) -> randindex(x, y)[3]
            elseif loss == "omARI" 
                lossfn = (x::ClustLabelVector, y::ClustLabelVector) -> 1 - randindex(x, y)[1]
            elseif loss == "VI"
                lossfn = varinfo
            elseif loss == "ID"
                lossfn = (x::ClustLabelVector, y::ClustLabelVector) -> infodist(x, y; normalised = false)
            end
        end
        clusts = samples.clusts
        lossmatrix = zeros(length(clusts), length(clusts))
        for i in eachindex(clusts)
            for j in (i+1):length(clusts)
                lossmatrix[i, j] = lossfn(clusts[i], clusts[j]) 
            end
        end
        lossmatrix += transpose(lossmatrix)
        i = argmin(vec(sum(lossmatrix, dims = 1)))
        return (clusts[i], i)
    end
end

"""
    binderloss(a::ClustLabelVector, b::ClustLabelVector;
    normalised = true) -> Float64

Computes the Binder loss between two clusterings. If `normalised = true` then the result is equal to one minus the rand index between `a` and `b`. 
"""
function binderloss(
    a::ClustLabelVector,
    b::ClustLabelVector;
    normalised = true
    )::Float64
    length(a) == length(b) || throw(ArgumentError("Length of the input vectors must be equal."))
    n = length(a)
    return randindex(a, b)[3] * (normalised ? 1 : binomial(n, 2))
end

@doc raw"""
    infodist(a::ClustLabelVector, b::ClustLabelVector;
    normalised = true) -> Float64

Computes the information distance between two clusterings. The information distance is defined as 
```math 
d_{\mathrm{ID}}(a, b) = \max \{H(A), H(B)\} - I(A, B)
```
where ``A`` and ``B`` are the cluster membership probability functions for ``a`` and ``b`` respectively, ``H`` denotes the entropy of a distribution, and ``I`` denotes the mutual information between two distributions. The information distance has range ``[0, \log N]`` where ``N`` is the number of observations. If normalised = true, the result is scaled to the range ``[0, 1]``.
"""
function infodist(
    a::ClustLabelVector,
    b::ClustLabelVector;
    normalised = true
    )::Float64
    length(a) != length(b) && throw(ArgumentError("Length of the input vectors must be equal."))
    n = length(a)
    hu = entropy(counts(a)/n)
    hv = entropy(counts(b)/n)
    return normalised ? 1 - mutualinfo(a, b; normed = false) / maximum([hu, hv]) : maximum([hu, hv]) - mutualinfo(a, b; normed = false)
end