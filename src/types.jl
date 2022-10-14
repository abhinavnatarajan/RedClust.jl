using Match, Printf, Dates

const ClustLabelType = Int
"Alias for Vector{Int}"
const ClustLabelVector = Vector{ClustLabelType}

const _numGibbs_default = 5
const _numMH_default = 1

"""
    MCMCOptionsList(; 
    numiters = 5000, 
    burnin = floor(0.2 * numiters), 
    thin = 1, 
    numGibbs = 5, 
    numMH = 1)

List of options for running the MCMC. 

# Constructor arguments
- `numiters::Integer`: number of iterations to run.
- `burnin::Integer`: number of iterations to discard as burn-in.
- `thin::Integer`: will keep every `thin` samples.
- `numGibbs:Integer`: number of intermediate Gibbs scans in the split-merge step.
- `numMH:Integer`: number of split-merge steps per MCMC iteration.
"""
struct MCMCOptionsList
    numiters::Int
    burnin::Int
    thin::Int
    numGibbs::Int
    numMH::Int
    numsamples::Int
    function MCMCOptionsList(;
        numiters::Int = 5000, 
        burnin::Int = Integer(floor(0.2 * numiters)),
        thin::Int = 1,
        numGibbs::Int = _numGibbs_default,
        numMH::Int = _numMH_default)
        # input validation
        if numiters < 1
            error("numiters must be ≥ 1.")
        end
        if burnin > numiters
            error("burnin must be < numiters")
        end
        if thin < 1
            error("thin must be positive.")
        end
        if numGibbs < 0
            error("numGibbs must be non-negative.")
        end
        if numMH < 0
            error("numMH must be non-negative.")
        end
        numsamples = Int(floor((numiters - burnin) / thin))
        new(numiters, burnin, thin, numGibbs, numMH, numsamples)
    end
end

function Base.show(io::IO, ::MIME"text/plain", options::MCMCOptionsList)
    printstyled(io, "MCMC Options\n"; color = :green, bold = true)
    println(io, "$(options.numiters) iteration$(options.numiters != 1 ? "s" : "")")
    println(io, "$(options.burnin) burnin iteration$(options.burnin != 1 ? "s" : "")")
    println(io, "$(options.numsamples) sample$(options.numsamples != 1 ? "s" : "")")
    println(io, "$(options.numGibbs) restricted Gibbs step$(options.numGibbs != 1 ? "s" : "") per split-merge step")
    print(io, "$(options.numMH) split-merge step$(options.numMH != 1 ? "s" : "") per iteration")
end

@doc raw"""
    PriorHyperparamsList(; [kwargs])

Contains the prior hyperparameters for the model. 

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
- `maxK::Int = 0`: maxmimum number of clusters to allow. If set to 0 then no maximum is imposed.

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
    proposalsd_r::Float64 = sqrt(η) / σ
	u::Float64 = 1
	v::Float64 = 1
    K_initial::Int = 1
    repulsion::Bool = true
    maxK::Int = 0
end

function Base.show(io::IO, ::MIME"text/plain", params::PriorHyperparamsList)
    printstyled(io, "Model Hyperparameters\n"; color = :green, bold = true)
    printstyled(io, "Likelihood Hyperparameters\n"; color = :magenta, bold = true)
    println(io, "δ₁ = $(prettynumber(params.δ1))")
    println(io, "δ₂ = $(prettynumber(params.δ2))")
    println(io, "α = $(prettynumber(params.α))")
    println(io, "β = $(prettynumber(params.β))")
    println(io, "ζ = $(prettynumber(params.ζ))")
    println(io, "γ = $(prettynumber(params.γ))")
    printstyled(io, "Partition Prior Hyperparameters\n"; color = :magenta, bold = true)
    println(io, "η = $(prettynumber(params.η))")
    println(io, "σ = $(prettynumber(params.σ))")
    println(io, "u = $(prettynumber(params.u))")
    println(io, "v = $(prettynumber(params.v))")
    printstyled(io, "Miscellaneous Hyperparameters\n"; color = :magenta, bold = true)
    println(io, "Proposal standard deviation for sampling r = $(prettynumber(params.proposalsd_r))")
    println(io, "Repulsion is enabled? $(params.repulsion)")
    println(io, "Maximum number of clusters = $(params.maxK > 0 ? params.maxK : "none")")
    print(io, "Initial number of clusters = $(params.K_initial)")
end

Base.@kwdef mutable struct MCMCState
    clusts::ClustLabelVector 
    r::Float64
    p::Float64
    clustsizes::Vector{Int} = counts(clusts, 1:length(clusts)) 
    K::Int = sum(clustsizes .> 0)
end

"""
    MCMCData(points::AbstractVector{<:AbstractVector{<:Float64}})    
    MCMCData(dissimilaritymatrix::AbstractMatrix{Float64})
    
Contains the pairwise dissimilarities for the MCMC sampler.  
"""
struct MCMCData 
    D::Matrix{Float64}
    logD::Matrix{Float64}
    function MCMCData(D::AbstractMatrix{Float64}) 
        if any(D .!= D')
            error("D must be symmetric.")
        end
        if size(D, 1) != size(D, 2)
            error("D must be a square matrix.")
        end
        new(Matrix(D), Matrix(log.(D .+ I(size(D,1)))))
    end
end

function MCMCData(pnts::AbstractVector{<:AbstractVector{<:Float64}})
    D = pairwise(Euclidean(), makematrix(pnts), dims = 2)
    MCMCData(D)
end

function Base.show(io::IO, data::MCMCData)
    printstyled(io, "MCMC data : "; color = :green, bold = true)
    print(io, "$(size(data.D, 1))×$(size(data.D, 1)) dissimilarity matrix.")
end

"""
Struct containing MCMC samples. 
# Fields
- `options::MCMCOptionsList`: options passed to the sampler. 
- `params::PriorHyperparamsList`: prior hyperparameters used by the sampler. 
- `clusts::Vector{ClustLabelVector}`: contains the clustering allocations. `clusts[i]` is a vector containing the clustering allocation from the `i`th sample. 
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
- `logposterior::Vector{Float64}`: a function proportional to the log-posterior for each sample, with constant of proportionality equal to the normalising constant of the partition prior.
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
    logposterior::Vector{Float64}
    options::MCMCOptionsList
    params::PriorHyperparamsList
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
        x.logposterior = zeros(numsamples)
        x.splitmerge_acceptances = Vector{Bool}(undef, numiters * numMH)
        x.splitmerge_splits = Vector{Bool}(undef, numiters * numMH)
        x.r_acceptances = Vector{Bool}(undef, numiters)
        return x
    end
end

function Base.show(io::IO, ::MIME"text/plain", result::MCMCResult)
    printstyled(io, "MCMC Summary\n"; color = :green, bold = true)
    printstyled(io, "General\n"; color = :magenta, bold = true)
    println(io, "$(result.options.numiters) iteration$(result.options.numiters != 1 ? "s" : "")")
    println(io, "$(result.options.burnin) iteration$(result.options.numiters != 1 ? "s" : "") discarded as burnin")
    println(io, "$(result.options.numsamples) sample$(result.options.numsamples != 1 ? "s" : "")")
    println(io, "$(result.options.numGibbs) restricted Gibbs step$(result.options.numGibbs != 1 ? "s" : "") per split-merge step")
    println(io, "$(result.options.numMH) split-merge step$(result.options.numMH != 1 ? "s" : "") per iteration")
    println(io, "Acceptance rate for split-merge steps = $(prettynumber(result.splitmerge_acceptance_rate))")
    println(io, "Acceptance rate for sampling r = $(prettynumber(result.r_acceptance_rate))")
    println(io, "Runtime = $(prettytime(result.runtime))")
    println(io, "Time per iteration : $(prettytime(result.mean_iter_time))")
    println(io, "")
    printstyled(io, "Summary for K\n"; color = :magenta, bold = true)
    println(io, "IAC : $(prettynumber(result.K_iac))")
    println(io, "ESS : $(prettynumber(result.K_ess))")
    println(io, "ESS per sample : $(prettynumber(result.K_ess / result.options.numsamples))")
    println(io, "Posterior mean : $(prettynumber(result.K_mean))")
    println(io, "Posterior variance : $(prettynumber(result.K_variance))")
    println(io, "")
    printstyled(io, "Summary for r\n"; color = :magenta, bold = true)
    println(io, "IAC : $(prettynumber(result.r_iac))")
    println(io, "ESS : $(prettynumber(result.r_ess))")
    println(io, "ESS per sample : $(prettynumber(result.r_ess / result.options.numsamples))")
    println(io, "Posterior mean : $(prettynumber(result.r_mean))")
    println(io, "Posterior variance : $(prettynumber(result.r_variance))")
    println(io, "")
    printstyled(io, "Summary for p\n"; color = :magenta, bold = true)
    println(io, "IAC : $(prettynumber(result.p_iac))")
    println(io, "ESS : $(prettynumber(result.p_ess))")
    println(io, "ESS per sample : $(prettynumber(result.p_ess / result.options.numsamples))")
    println(io, "Posterior mean : $(prettynumber(result.p_mean))")
    print(io, "Posterior variance : $(prettynumber(result.p_variance))")
end

function Base.show(io::IO, result::MCMCResult)
    print(io, "MCMC result with $(result.options.numsamples) samples")
end

function prettytime(t)
    if t < 1e-6
        x = Dates.canonicalize(Dates.Nanosecond(Integer(floor(t * 1e9))))
    elseif t ≥ 1e-6 && t < 1e-3
        x = Dates.canonicalize(Dates.Microsecond(Integer(floor(t * 1e6))))
    elseif t ≥ 1e-3 && t < 1
        x = Dates.canonicalize(Dates.Millisecond(Integer(floor(t * 1e3))))
    else
        x = Dates.canonicalize(Dates.Second(Integer(floor(t)))) 
    end
    out = @match x.periods[1] begin
        _::Nanosecond => @sprintf("%.3f ns", t * 1e9)
        _::Microsecond => @sprintf("%.3f μs", t * 1e6)
        _::Millisecond => @sprintf("%.3f μs", t * 1e3)
        _::Second => @sprintf("%.3f s", t)
        _::Minute => @sprintf("%dm%ds", x.periods[1].value, x.periods[2].value)
        _::Hour => @sprintf("%dh%dm%ds", x.periods[1].value, x.periods[2].value, x.periods[3].value)
        _::Day => @sprintf("%dD%dh%dm", x.periods[1].value, x.periods[2].value, x.periods[3].value)
    end
    out
end

prettynumber(x) = x < 1 ? @sprintf("%.3e", x) : @sprintf("%.3f", x)
