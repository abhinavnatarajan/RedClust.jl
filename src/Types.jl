const ClustLabelType = Int
"Alias for Vector{Int}"
const ClustLabelVector = Vector{ClustLabelType}

const _numGibbs_default = 5
const _numMH_default = 1

"""
List of options for running the MCMC. 

Constructor:
    MCMCOptionsList(; [numiters, burnin, thin, numGibbs, numMH, pointestimation, usesalso])

# Constructor arguments
- `numiters::Integer = 5000`: number of iterations to run.
- `burnin::Integer = floor(0.2 * numiters)`: number of iterations to discard as burn-in.
- `thin::Integer = 1`: will keep every `thin` samples.
- `numGibbs:Integer = 1`: number of intermediate Gibbs scans in the split-merge step.
- `numMH:Integer = 1`: number of split-merge steps per MCMC iteration.
- `pointestimation::Bool = true`: whether to compute a point-estimate clustering from the posterior samples.
- `usesalso::Bool = false`: whether to use the SALSO algorithm to compute a point estimate.
"""
Base.@kwdef struct MCMCOptionsList
    numiters::Int = 5000
    burnin::Int = Integer(floor(0.2 * numiters))
    thin::Int = 1
    numGibbs::Int = _numGibbs_default
    numMH::Int = _numMH_default
    pointestimation::Bool = true
    usesalso::Bool = false 
    numsamples::Int = Int(floor((numiters - burnin) / thin))
end

@doc raw"""
Contains the prior hyperparameters for the model. Call the constructor with the line 
    
    PriorHyperparamsList(; [kwargs])

# Constructor arguments
- `δ1::Float64 = 1`: the parameter ``\delta_1``.
- `α::Float64 = 1`: the shape parameter ``\alpha`` in the prior for each ``\lambda_k``.
- `β::Float64 = 1`: the rate parameter ``\beta`` in the prior for each ``\lambda_k``.
- `δ2::Float64 = 1`: the parameter ``\delta_2``.
- `ζ::Float64 = 1`: the shape parameter ``\zeta`` in the prior for each ``\theta_{kt}``.
- `γ::Float64 = 1`: the rate parameter ``\gamma`` in the prior for each ``\theta_{kt}``.
- `η::Float64 = 1`: the shape parameter ``\eta`` in the prior for ``r``.
- `σ::Float64 = 1`: the rate parameter ``\sigma`` in the prior for ``r``.
- `u::Float64 = 1`: the parameter ``u`` in the prior for ``p``.
- `v::Float64 = 1`: the parameter ``v`` in the prior for ``p``.
- `proposalsd_r::Float64 = 1`: standard deviation of the truncated Normal proposal when sampling `r` via a Metropolis-Hastings step.
- `K_initial::Int = 1`: initial number of clusters to start the MCMC.
- `repulsion::Bool = true`: whether to use the repulsive component of the likelihood. 
- `maxK::Int = 0`: maxmimum number of clusters to allow. If set to 0 then we can get up to `n` clusters, where `n` is the number of observations.

# See also
[`fitprior`](@ref)
"""
Base.@kwdef mutable struct PriorHyperparamsList
	δ1::Float64 = 1
    δ2::Float64 = 1
	α::Float64 = 1
	β::Float64 = 1
	ζ::Float64 = 1
	γ::Float64 = 1
	η::Float64 = 1
	σ::Float64 = 1
    proposalsd_r::Float64 = 1
	u::Float64 = 1
	v::Float64 = 1
    K_initial::Int = 1
    repulsion::Bool = true
    maxK::Int = 0
end

Base.@kwdef mutable struct MCMCState
    clusts::ClustLabelVector 
    r::Float64
    p::Float64
    clustsizes::Vector{Int} = counts(clusts, 1:length(clusts)) 
    K::Int = sum(clustsizes .> 0)
end

"Contains the pairwise dissimilarities. Construct by calling `MCMCData(D)`, where `D` is a square matrix of pairwise dissimilarities."
Base.@kwdef struct MCMCData 
    D::Matrix{Float64}
    logD::Matrix{Float64} = log.(D + I(size(D,1)))
end

"""
Struct containing MCMC samples. 
# Fields
- `options::MCMCOptionsList`: options passed to the sampler. 
- `params::PriorHyperparamsList`: prior hyperparameters used by the sampler. 
- `clusts::Vector{ClustLabelVector}`: contains the clustering allocations. `clusts[i]` is a vector containing the clustering allocation from the `i`th sample. 
- `pntestimate::ClustLabelVector`: contains the point-estimate clustering allocation. If `options.pointestimation == false` then this is a vector of zeros. 
- `posterior_coclustering::Matrix{Float64}`: the posterior coclustering matrix. 
- `K::Vector{Int}`: posterior samples of ``K``, i.e., the number of clusters. `K[i]` is the number of clusters in `clusts[i]`.
- `r::Vector{Float64}`: posterior samples of the parameter ``r``.
- `p::Vector{Float64}`: posterior samples of the parameter ``p``.
- `K_mean`, `r_mean`, `p_mean`: posterior mean of `K`, `r`, and `p` respectively. 
- `K_variance`, `r_variance`, `p_variance`: posterior variance of `K`, `r`, and `p` respectively. 
- `K_acf::Vector{Float64}`, `r_acf::Vector{Float64}`, `p_acf::Vector{Float64}`: autocorrelation function for `K`, `r`, and `p` respectively. 
- `K_iac`, `r_iac`, and `p_iac`: integrated autocorrelation coefficient for `K`, `r`, and `p` respectively. 
- `K_ess::Float64`, `r_ess::Float64`, and `p_ess::Float64`: effective sample size for `K`, `r`, and `p` respectively.
- `loglik::Vector{Float64}`: log-likelihood for each sample. 
- `splitmerge_splits`: Boolean vector indicating the iterations when a split proposal was used in the split-merge step. Has length `numMH * numiters` (see [`MCMCOptionsList`](@ref)).
- `splitmerge_acceptance_rate`: acceptance rate of the split-merge proposals. 
- `r_acceptances`: Boolean vector indicating the iterations (including burnin and the thinned out iterations) where the Metropolis-Hastings proposal for `r` was accepted. 
- `r_acceptance_rate:`: Metropolis-Hastings acceptance rate for `r`.
- `runtime`: total runtime for all iterations.
- `mean_iter_time`: average time taken for each iteration. 
"""
Base.@kwdef mutable struct MCMCResult
    clusts::Vector{ClustLabelVector}
    posterior_coclustering::Matrix{Float64}
    K::Vector{Int}
    K_ess::Float64
    K_acf::Vector{Float64}
    K_iac::Float64
    K_mean::Float64
    K_variance::Float64
    r::Vector{Float64}
    r_ess::Float64
    r_acf::Vector{Float64}
    r_iac::Float64
    r_mean::Float64
    r_variance::Float64
    p::Vector{Float64}
    p_ess::Float64
    p_acf::Vector{Float64}
    p_iac::Float64
    p_mean::Float64
    p_variance::Float64
    splitmerge_acceptances::Vector{Bool}
    r_acceptances::Vector{Bool}
    r_acceptance_rate::Float64
    splitmerge_splits::Vector{Bool}
    splitmerge_acceptance_rate::Float64
    runtime::Float64
    mean_iter_time::Float64
    loglik::Vector{Float64}
    options::MCMCOptionsList
    params::PriorHyperparamsList
    pntestimate::Vector{Int}
    function MCMCResult(data::MCMCData, options::MCMCOptionsList, params::PriorHyperparamsList)
        x = new()
        numpts = size(data.D, 1)
        x.options = options
        x.params = params
        numiters = options.numiters
        numsamples = options.numsamples
        numMH = options.numMH
        x.clusts = [zeros(ClustLabelType, numpts) for i in 1:numsamples]
        x.posterior_coclustering = zeros(numpts, numpts)
        x.K = Int.(zeros(numsamples))
        x.K_acf = zeros(numsamples)
        x.r = zeros(numsamples)
        x.r_acf = zeros(numsamples)
        x.p = zeros(numsamples)
        x.p_acf = zeros(numsamples)
        x.loglik = zeros(numsamples)
        x.splitmerge_acceptances = Vector{Bool}(undef, numiters * numMH)
        x.splitmerge_splits = Vector{Bool}(undef, numiters * numMH)
        x.r_acceptances = Vector{Bool}(undef, numiters)
        x.pntestimate = zeros(Int, numpts)
        return x
    end
end