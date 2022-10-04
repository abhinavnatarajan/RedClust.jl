using Random: rand
using StatsBase: sample, counts, mean_and_var, mean, var
using Distributions: Normal, Uniform, Gamma, Beta, logpdf, truncated
using SpecialFunctions: loggamma
using LinearAlgebra: I
using Clustering: kmedoids, varinfo
using ProgressBars: ProgressBar
using RCall: rcopy, @R_str



function loglik(
    data::MCMCData, 
    state::MCMCState,
    params::PriorHyperparamsList,
    clustering::Bool,
    partition::Bool
    )::Float64
    D = data.D
    logD = data.logD
    clusts = state.clusts
    clustsizes = state.clustsizes
    K = state.K
    r = state.r
    p = state.p
    δ1 = params.δ1
    δ2 = params.δ2
    α = params.α
    β = params.β
    ζ = params.ζ
    γ = params.γ
    η = params.η
    σ = params.σ
    u = params.u
    v = params.v
    repulsion = params.repulsion

	C = findall(clustsizes .> 0)
    n = size(D, 1)
	L = 0
	
	if clustering 
		# Cohesive part of likelihood 
		for k = 1:K
			clust_k = clusts .== C[k]
			sz_k = clustsizes[C[k]]
			pairs_k = numpairs(sz_k)
			a = α + δ1 * pairs_k
			b = β + sum(D[clust_k, clust_k]) / 2 
			L += (δ1 - 1) * sum(logD[clust_k, clust_k]) / 2 - pairs_k * loggamma(δ1) +
				α * log(β) - loggamma(α) +
				loggamma(a)  - a * log(b)
        end
		
		# Repulsive part of likelihood
		if repulsion
			for k = 1:K
				clust_k = clusts .== C[k]
				sz_k = clustsizes[C[k]]
				for t = (k + 1):K
					clust_t = clusts .== C[t]
					sz_t = clustsizes[C[t]]	
					pairs_kt = sz_k * sz_t
					z = ζ + δ2 * pairs_kt
					g = γ + sum(D[clust_k, clust_t])
					L += (δ2 - 1) * sum(logD[clust_k, clust_t]) - pairs_kt * loggamma(δ2) +
						ζ * log(γ) - loggamma(ζ) +
						loggamma(z) - z * log(g)
                end
            end
        end
    end
	
	if partition
        L += loggamma(K + 1) + (n - K) * log(p) + (r * K) * log(1 - p) - K * loggamma(r)
        for nj in clustsizes[C]
            L += log(nj) + loggamma(nj + r - 1)
        end
    end
	return L
end

function sample_r!(state::MCMCState, params::PriorHyperparamsList)
    r = state.r
    p = state.p
    clustsizes = state.clustsizes
    C = clustsizes[findall(clustsizes .> 0)]
    K = state.K
    η = params.η
    σ = params.σ
    proposalsd_r = params.proposalsd_r
    temp = sample_r(r, p, C, K, η, σ, proposalsd_r)
    state.r = temp.r
    return (r = temp.r, accept = temp.accept)
end

function sample_r(
    r::Float64, 
    p::Float64, 
    C::Vector{Int},
    K::Integer, 
    η::Float64,
    σ::Float64,
    proposalsd_r::Float64
    )
	# proposal
    proposal_distribution = truncated(
        Normal(r, proposalsd_r), 
        lower = 0,
        upper = Inf
    )
    r_candidate = rand(proposal_distribution)
    reverse_distribution = truncated(
        Normal(r_candidate, proposalsd_r), 
        lower = 0,
        upper = Inf
    )
    
	# calculate transition probability
	logprobcandidate = (η-1) * log(r_candidate) + K * (r_candidate * log(1-p) - loggamma(r_candidate)) - r_candidate * σ
	logprobcurrent = (η-1) * log(r) + K * (r * log(1-p) - loggamma(r)) - r * σ 
	for nk in C
		logprobcandidate = logprobcandidate + loggamma(nk - 1 + r_candidate)
		logprobcurrent = logprobcurrent + loggamma(nk - 1 + r)
    end
	# Log-ratio of proposal densities
	logproposalratio = logpdf(proposal_distribution, r_candidate) -
    logpdf(reverse_distribution, r)
		
	rnew = r
	accept = false
	# Metropolis acceptance step
	if log(rand(Uniform())) < 
        minimum([0, logprobcandidate - logprobcurrent - logproposalratio]) # accept
		rnew = r_candidate
		accept = true
    end
	return (r = rnew, accept = accept)
end

function sample_p!(state::MCMCState, params::PriorHyperparamsList)::Float64
    K = state.K
    n = length(state.clusts)
    r = state.r
    u = params.u
    v = params.v
    state.p = sample_p(K, n, r, u, v)
end

function sample_p(
    K::Int, 
    n::Int,
    r::Float64, 
    u::Float64, 
    v::Float64
    )::Float64
    rand(Beta(n - K + u, r * K + v))
end

## Sample clustering allocation labels via Gibbs sampling
function sample_labels_Gibbs!(
    data::MCMCData,
    state::MCMCState,
    params::PriorHyperparamsList,
    items_to_reallocate::Vector{Int} = Vector(1:length(state.clusts)),
    candidate_clusts::Vector{ClustLabelType} = ClustLabelType[],
    final_clusts::ClustLabelVector = ClustLabelVector()
    )

	D = data.D
    logD = data.logD
    clusts = state.clusts
    clustsizes = state.clustsizes
	n = length(clusts)
    r = state.r
    p = state.p
    δ1 = params.δ1
    δ2 = params.δ2
    α = params.α
    β = params.β
    ζ = params.ζ
    γ = params.γ
    repulsion = params.repulsion
    maxK = params.maxK
	free_label_set = isempty(candidate_clusts) # free allocation of cluster labels is allowed
	free_allocation = isempty(final_clusts) # force the result of the Gibbs scan (when we want to compute transition probability only)
	log_transition_prob = 0 # transition probability of the Gibbs scan

	α_i = zeros(n)
    β_i = zeros(n)
    ζ_i = zeros(n)
    γ_i = zeros(n)
    sum_logD_i = zeros(n)
    L2_ik_prime = zeros(n)
    αβratio = α * log(β) - loggamma(α)
    ζγratio = ζ * log(γ) - loggamma(ζ)
	for i in items_to_reallocate
		clustsizes[clusts[i]] = clustsizes[clusts[i]] - 1
		clusts[i] = -1
		C_i = findall(clustsizes .> 0)
		K_i = length(C_i)

		candidate_clusts = C_i
        if free_label_set && (maxK == 0 || K_i < maxK) && K_i < n
			candidate_clusts = append!(candidate_clusts, findall(clustsizes .== 0)[1])
        end
		m = length(candidate_clusts)
		
		# Calculate everything that we will need
		
		for k in C_i
			clust_k = (clusts .== k)
			sz_clust_k = clustsizes[k]
			α_i[k] = α + δ1 * sz_clust_k
			β_i[k] = β + sum(D[i, clust_k])
			ζ_i[k] = ζ + δ2 * sz_clust_k
			γ_i[k] = γ + sum(D[i, clust_k])
			sum_logD_i[k] = sum(logD[i, clust_k])
        end
		
		L1 = zeros(m)
		L2 = zeros(m)
		log_prior_ratio = zeros(m)
		logprobs = zeros(m)
		
		for k in 1:m # compute cohesive likelihood component for each proposed cluster
			if clustsizes[candidate_clusts[k]] == 0 # new cluster
				log_prior_ratio[k] = log(K_i + 1) + r * log(1 - p)
				L1[k] = 0
			else 
				sz_candidate_clust_k = clustsizes[candidate_clusts[k]]
				L1[k] = loggamma(α_i[candidate_clusts[k]]) + αβratio - 
                α_i[candidate_clusts[k]] * log(β_i[candidate_clusts[k]]) +
					(δ1 - 1) * sum_logD_i[candidate_clusts[k]] - sz_candidate_clust_k * loggamma(δ1)
                log_prior_ratio[k] = log(sz_candidate_clust_k + 1) + log(p) + log(sz_candidate_clust_k - 1 + r) - log(sz_candidate_clust_k)
            end
		end
		
		logprobs .= log_prior_ratio .+ L1
		
		if repulsion
            for t in C_i 
                L2_ik_prime[t] = loggamma(ζ_i[t]) - ζ_i[t] * log(γ_i[t]) +
                ζγratio + (δ2 - 1) * sum_logD_i[t] - clustsizes[t] * loggamma(δ2)
            end
            L2_i = sum(L2_ik_prime[C_i])
            for k = 1:m
                if clustsizes[candidate_clusts[k]] == 0
                    L2[k] = L2_i
                else
                    L2[k] = L2_i - L2_ik_prime[candidate_clusts[k]]
                end
            end
			logprobs .+= L2
		end 
		if free_allocation 
			k = sample_logweights(logprobs)
			ci_new = candidate_clusts[k]
        else 
			ci_new = final_clusts[i]
			k = findall(candidate_clusts .== ci_new)[1]
		end
		
		clusts[i] = ci_new
		clustsizes[ci_new] = clustsizes[ci_new] + 1
		
		# calculate the transition probability
		logprobs .+= minimum(logprobs)
		probs = exp.(logprobs)
		probs ./= sum(probs)
		log_transition_prob += log(probs[k])
    end
    state.K = sum(clustsizes .> 0)
    return log_transition_prob
end

function sample_labels!(
    data::MCMCData,
    state::MCMCState,
    params::PriorHyperparamsList,
    options::MCMCOptionsList
    ) 
    # Unpack
    
    n = length(state.clusts) 
    r = state.r
    p = state.p
    numMH = options.numMH
    numGibbs = options.numGibbs

	accept = fill(false, numMH)
    split = fill(false, numMH)
    for MH_counter in 1:numMH
        clusts = state.clusts
        clustsizes = state.clustsizes
        K = state.K
        # BEGIN MH
        # Pick chaperones
        i, j = sample(1:n, 2, replace = false)
        ci = clusts[i] # labels of their current cluster allocations
        cj = clusts[j]
        
        # automatically reject if we have max number of clusters
        if params.maxK > 0 && ci == cj && length(findall(clustsizes .> 0)) >= params.maxK
            continue
        end

        # choose elements to reallocate
        S = findall(clusts .== ci .|| clusts .== cj)
        S = S[S .!= i .&& S .!= j]
        
        # initialise random launch state
        claunch = copy(clusts)
        szlaunch = copy(clustsizes)
        Klaunch = K
        if ci == cj 
            claunch[i] = findall(clustsizes .== 0)[1]
            szlaunch[ci] = szlaunch[ci] - 1
            szlaunch[claunch[i]] = szlaunch[claunch[i]] + 1
            Klaunch = K + 1
        end
        candidate_clusts = [claunch[i], claunch[j]]
        for k in S 
            claunch[k] = sample(candidate_clusts)
            szlaunch[clusts[k]] = szlaunch[clusts[k]] - 1
            szlaunch[claunch[k]] = szlaunch[claunch[k]] + 1
        end
        launchstate = MCMCState(claunch, r, p, szlaunch, Klaunch)
        
        # restricted Gibbs scans
        for scan_counter in 1:numGibbs 
            sample_labels_Gibbs!(data, 
            launchstate, params, S, candidate_clusts) 
        end
        
        if ci == cj # split proposal
            split[MH_counter] = true
            # final Gibbs scan
            log_transition_prob = sample_labels_Gibbs!(data, 
            launchstate, params, S, candidate_clusts) 
            finalstate = launchstate
            cfinal = finalstate.clusts
            szfinal = finalstate.clustsizes
            Kfinal = finalstate.K
            
            # log prior ratio for MH acceptance ratio
            log_prior_ratio = log(K + 1) + r * log(1 - p) - log(p) - loggamma(r) + 
                loggamma(szfinal[cfinal[i]] - 1 + r) + loggamma(szfinal[cfinal[j]] - 1 + r) +
                log(szfinal[cfinal[i]]) + log(szfinal[cfinal[j]]) + 
                -(loggamma(clustsizes[ci] - 1 + r) + log(clustsizes[ci])) 
            
            # proposal density ratio pr(new | old) / pr(old | new)
            log_proposal_ratio = log_transition_prob
            # note that the reverse proposal density is just 1
        else # merge
            cfinal = copy(launchstate.clusts)
            szfinal = copy(launchstate.clustsizes)
            Kfinal = launchstate.K
            clust_i = findall(cfinal .== ci)
            sz_clust_i = length(clust_i)
            cfinal[clust_i] .= cj
            szfinal[ci] = 0
            szfinal[cj] = szfinal[cj] + sz_clust_i
            Kfinal -= 1
            finalstate = MCMCState(cfinal, r, p, szfinal, Kfinal)

            # log prior ratio for MH acceptance ratio
            log_prior_ratio = -(log(K) + r * log(1 - p) - log(p) - loggamma(r)) + 
            loggamma(szfinal[cj] - 1 + r) + log(szfinal[cj]) +
            -(loggamma(clustsizes[ci] - 1 + r) + loggamma(clustsizes[cj] - 1 + r)
            + log(clustsizes[ci]) + log(clustsizes[cj]))
            
            # proposal density ratio pr(new | old) / pr(old | new)
            log_transition_prob = sample_labels_Gibbs!(data, 
            launchstate, params, S, candidate_clusts, clusts)
                                
            log_proposal_ratio = -log_transition_prob
            # note that the reverse proposal density is just 1
        end
        
        # likelihood ratio for MH acceptance ratio
        log_lik_ratio = loglik(data,
            finalstate, params, true, false) - 
            loglik(data, state, params, true, false)
        
        # MH acceptance step
        log_acceptance_ratio = minimum(
            [0, log_prior_ratio + log_lik_ratio - log_proposal_ratio])
        if (log(rand(Uniform())) < log_acceptance_ratio) # accept
            state = finalstate
            accept[MH_counter] = true
        end
        # END MH
    end
    
	
	# Final Gibbs scan
    temp = sample_labels_Gibbs!(data, state, params)
    return (accept = accept, split = split)
end

"""
    runsampler(data[, options, params])

Runs the MCMC sampler on the data.

# Arguments
- `data::MCMCData`: contains the distance matrix. 
- `options = MCMCOptionsList()`: contains the number of iterations, burnin, etc. 
- `params = PriorHyperparamsList()`: contains the prior hyperparameters for the model.

# Returns
A struct of type `MCMCResult` containing the MCMC samples, convergence diagnostics, and summary statistics.

# See also
[`MCMCData`](@ref), [`MCMCOptionsList`](@ref), [`fitprior`](@ref), [`PriorHyperparamsList`](@ref).
"""
function runsampler(data::MCMCData, 
    options::MCMCOptionsList = MCMCOptionsList(), 
    params::Union{PriorHyperparamsList, Nothing} = nothing,
    init::Union{MCMCState, Nothing} = nothing
    )::MCMCResult
    numiters = options.numiters
    burnin = options.burnin
    thin = options.thin
    numsamples = options.numsamples
    if isnothing(params)
        params = fitprior(data.D, "k-medoids", true).params
    end
    if isnothing(init)
        init = MCMCState(
        clusts = kmedoids(data.D, 
            (params.maxK > 0 ? minimum([params.maxK, params.K_initial]) : params.K_initial); 
            maxiter=1000).assignments,
        r = rand(Gamma(params.η, 1 / params.σ)),
        p = rand(Beta(params.u, params.v))
        )
    end
    result = MCMCResult(data, options, params)
    state = init

    # Start sampling
    j::Int = 1
    print("Start sampling...\n")
    runtime = 0
    runtime = @elapsed (
    for i in ProgressBar(1:numiters)
        result.r_acceptances[i] = sample_r!(state, params).accept
        sample_p!(state, params)
        temp = sample_labels!(data, state, params, options)
        temprange = ((i-1) * options.numMH + 1) : (i * options.numMH)
        result.splitmerge_acceptances[temprange] .= temp.accept
        result.splitmerge_splits[temprange] .= temp.split

        # Record sample if necessary
        if i > burnin && (i - burnin) % thin == 0
            result.clusts[j] .= sortlabels(state.clusts)
            result.K[j] = state.K
            result.r[j] = state.r
            result.p[j] = state.p
            result.loglik[j] = loglik(data, state, params, true, true)
            j = j + 1
        end
    end
    )

    print("Computing summary statistics and diagnostics...\n")
    # Create co-clustering matrix and compute its IAC
    result.posterior_coclustering = sum(adjacencymatrix.(result.clusts)) ./ numsamples

    # Get point estimate
    if options.pointestimation
        result.pntestimate = pointestimate(result.clusts; loss = "VI", usesalso = options.usesalso)
    end


    # Compute normalised autocorrelation and integrated autocorrelation coefficient of the samples
    # For K
    result.K_iac, result.K_ess, result.K_acf = iac_ess_acf(result.K)
    result.K_mean, result.K_variance = mean_and_var(result.K)

    # For r
    result.r_iac, result.r_ess, result.r_acf = iac_ess_acf(result.r)
    result.r_mean, result.r_variance = mean_and_var(result.r)

    # For p
    result.p_iac, result.p_ess, result.p_acf = iac_ess_acf(result.p)
    result.p_mean, result.p_variance = mean_and_var(result.p)

    # Acceptance rate of the MH steps
    if options.numMH > 0
        result.splitmerge_acceptance_rate = mean(result.splitmerge_acceptances)
    else
        result.splitmerge_acceptance_rate = 0
    end
    result.r_acceptance_rate = mean(result.r_acceptances)

    # Miscellaneous
    result.options = options
    result.params = params
    result.runtime = runtime
    result.mean_iter_time = runtime / numiters
    return result
end

function sample_rp(
    clustsizes::Vector{Int},
    options::MCMCOptionsList = MCMCOptionsList(),
    params::PriorHyperparamsList = PriorHyperparamsList()
    ) 
    # Unpack parameters
    numiters = options.numiters
    burnin = options.burnin
    thin = options.thin
	η = params.η
    σ = params.σ
    proposalsd_r = params.proposalsd_r
    u = params.u
    v = params.v
    C = clustsizes[findall(clustsizes .> 0)]
    n = sum(C)
    K = length(C)
	
	# Initialise result
    r = rand(Gamma(η, σ))
    p = rand(Beta(u, v))
	numsamples = Integer(floor((numiters - burnin) / thin))
	result = (r = zeros(numsamples), p = zeros(numsamples))
	j = 1
	for i in ProgressBar(1:numiters)
		# Sample r
		r = sample_r(r, p, C, K, η, σ, proposalsd_r).r
		# Sample p
		p = sample_p(K, n, r, u, v)
		
		# Record sample if necessary
		if i > burnin && (i - burnin) % thin == 0
			result.r[j] = r
			result.p[j] = p
			j = j + 1
        end
    end
	return result
end