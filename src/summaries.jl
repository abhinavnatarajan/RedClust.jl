using Clustering: randindex, varinfo, vmeasure, mutualinfo

@doc raw"""
    evaluateclustering(clusts::ClustLabelVector, truth::ClustLabelVector)
Returns a tuple of values quantifying the accuracy of the clustering assignments in `clusts` with respect to the ground truth clustering assignments in `truth`. 
# Return values
- `nbloss`: Normalised Binder loss (= 1 - rand index). Lower is better. 
- `ari`: Adjusted Rand index. Higher is better. 
- `vi`: Variation of Information distances. Lower is better. Range is ``[0, \log N]``.
- `nvi`: Normalised VI distance with range ``[0, 1]``.
- `id`: Information Distance. Lower is better. Range is ``[0, \log N]``.
- `nid`: Normalised Information Distance with range ``[0, 1]``.
- `vmeas`: V-Measure. Higher is better. 
- `nmi`: Normalised Mutual Information. Higher is better. 
"""
function evaluateclustering(clusts::ClustLabelVector, truth::ClustLabelVector)
    length(clusts) != length(truth) && throw(ArgumentError("Length of inputs must be equal."))
    n = length(clusts)
    ari, _, nbloss, _ = randindex(clusts, truth)
    vi = varinfo(clusts, truth)
    nvi = vi / log(n)
    id = infodist(clusts, truth; normalised = false)
    nid = id / log(n)
    vmeas = vmeasure(clusts, truth)
    nmi = mutualinfo(clusts, truth)
    return (nbloss = nbloss, ari = ari, vi = vi, nvi, id = id, nid = nid, vmeas = vmeas, nmi = nmi)
end

"""
    summarise(clusts::ClustLabelVector, truth::ClustLabelVector, printoutput = true)
Returns a string that summarises the clustering accuracy of `clusts` with respect to the ground truth in `truth`. If `printoutput = true`, the output will be printed to console before being returned. 
"""
function summarise(clusts::ClustLabelVector, truth::ClustLabelVector, printoutput=true)::String
    output = ""
    temp = evaluateclustering(clusts, truth)
    output *= "Clustering summary\n"
    output *= "Number of clusters : $(length(unique(clusts)))\n"
    output *= "Normalised Binder loss : $(temp.nbloss)\n"
    output *= "Adjusted Rand Index : $(temp.ari)\n"
    output *= "Normalised Variation of Information (NVI) distance : $(temp.nvi)\n"
    output *= "Normalised Information Distance (NID) : $(temp.nid)\n"
    output *= "Normalised Mutual Information : $(temp.nmi)\n"
    output *= "V-Measure : $(temp.vmeas)\n"
    if printoutput print(output) end
    return output
end

"""
    summarise(result::MCMCResult, printoutput = true)
Returns a string that summarises the MCMC output. If `printoutput = true`, the output will be printed to console before being returned. 
"""
function summarise(result::MCMCResult, printoutput=true)::String
    output = ""
    output *= "General Summary\n"
    output *= "Iterations : $(result.options.numiters)\n"
    output *= "Burn-in : $(result.options.burnin)\n"
    output *= "Thin : $(result.options.thin)\n"
    output *= "Number of samples : $(result.options.numsamples)\n"
    output *= "Number of intermediate Gibbs scans in each split-merge step : $(result.options.numGibbs)\n"
    output *= "Number of split-merge steps per iteration : $(result.options.numMH)\n"
    output *= "Acceptance rate for split-merge steps : $(result.splitmerge_acceptance_rate)\n"
    output *= "Acceptance rate for sampling r : $(result.r_acceptance_rate)\n"
    output *= "Runtime : $(result.runtime)\n"
    output *= "Time per iteration : $(result.mean_iter_time)\n"
    output *= "\n"
    output *= "Summary for K\n"
    output *= "IAC : $(result.K_iac)\n"
    output *= "ESS : $(result.K_ess)\n"
    output *= "ESS per sample : $(result.K_ess / result.options.numsamples)\n"
    output *= "Posterior mean : $(result.K_mean)\n"
    output *= "Posterior variance : $(result.K_variance)\n"
    output *= "\n"
    output *= "Summary for r\n"
    output *= "IAC : $(result.r_iac)\n"
    output *= "ESS : $(result.r_ess)\n"
    output *= "ESS per sample : $(result.r_ess / result.options.numsamples)\n"
    output *= "Posterior mean : $(result.r_mean)\n"
    output *= "Posterior variance : $(result.r_variance)\n"
    output *= "\n"
    output *= "Summary for p\n"
    output *= "IAC : $(result.p_iac)\n"
    output *= "ESS : $(result.p_ess)\n"
    output *= "ESS per sample : $(result.p_ess / result.options.numsamples)\n"
    output *= "Posterior mean : $(result.p_mean)\n"
    output *= "Posterior variance : $(result.p_variance)\n"
    if printoutput print(output) end
    return output
end