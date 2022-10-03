using Clustering: randindex, varinfo, vmeasure, mutualinfo

"""
    evaluateclustering(clusts::ClustLabelVector, truth::ClustLabelVector)
Returns a tuple of values quantifying the accuracy of the clustering assignments in `clusts` with respect to the ground truth clustering assignments in `truth`. 
# Return values
- `bloss`: Normalised Binder loss (= 1 - rand index). Lower is better. 
- `ari`: Adjusted Rand index. Higher is better. 
- `vidist`: Variation of Information distances. Lower is better. 
- `vmeas`: V-Measure. Higher is better. 
- `nmi`: Normalised Mutual Information. Higher is better. 
"""
function evaluateclustering(clusts::ClustLabelVector, truth::ClustLabelVector)
    bloss = binderloss(clusts, truth)
    ari = randindex(clusts, truth)[1]
    vidist = varinfo(clusts, truth)
    vmeas = vmeasure(clusts, truth)
    nmi = mutualinfo(clusts, truth)
    return (bloss, ari = ari, vidist = vidist, vmeas = vmeas, nmi = nmi)
end

function evaluateclustering(result::MCMCResult, truth::ClustLabelVector)
    reftruth = Ref(truth)
    bloss = binderloss.(result.clusts, reftruth)
    ari = (x -> randindex(x, truth)).(result.clusts, reftruth)
    vidist = varinfo.(result.clusts, reftruth)
    vmeas = vmeasure.(result.clusts, reftruth)
    nmi = mutualinfo.(result.clusts, reftruth)
    return (bloss, ari = ari, vidist = vidist, vmeas = vmeas, nmi = nmi)
end

function summarise(clusts::ClustLabelVector, truth::ClustLabelVector)::String
    output = ""
    temp = evaluateclustering(clusts, truth)
    output *= "Clustering summary\n"
    output *= "Number of clusters : $(length(unique(clusts)))\n"
    output *= "Normalised Binder loss : $(temp.bloss)\n"
    output *= "Variation of information distance : $(temp.vidist)\n"
    output *= "Adjusted rand index : $(temp.ari)\n"
    output *= "Normalised mutual information : $(temp.nmi)\n"
    output *= "V-Measure : $(temp.vmeas)\n"
    return output
end

"""
    summarise(result::MCMCResult, truth::ClustLabelVector = [], printoutput = true)
Returns a string that summarises the MCMC output. If `truth` is a vector of cluster labels, the output also contains clustering metrics to quantify the accuracy of the point-estimate from the MCMC with respect to `truth`. If `printoutput = true`, the output will be printed to console before being returned. 
"""
function summarise(result::MCMCResult, truth::ClustLabelVector=[], printoutput=true)::String
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
    if length(truth) > 0
        temp = evaluateclustering(result.pntestimate, truth)
        output *= "\n"
        output *= "Clustering summary\n"
        output *= "Number of clusters : $(length(unique(result.pntestimate)))\n"
        output *= "Normalised Binder loss : $(temp.bloss)\n"
        output *= "Variation of information distance : $(temp.vidist)\n"
        output *= "Adjusted rand index : $(temp.ari)\n"
        output *= "Normalised mutual information : $(temp.nmi)\n"
        output *= "V-Measure : $(temp.vmeas)\n"
    end
    if printoutput print(output) end
    return output
end

function binderloss(
    a::ClustLabelVector,
    b::ClustLabelVector
    )::Float64
    n = length(a)
    sum(abs.(adjacencymatrix(a) .- adjacencymatrix(b))) / (n * (n - 1))
end