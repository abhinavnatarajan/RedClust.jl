using Clustering: randindex, varinfo, vmeasure, mutualinfo

@doc raw"""
    evaluateclustering(clusts::ClustLabelVector, truth::ClustLabelVector)
Returns a named tuple of values quantifying the accuracy of the clustering assignments in `clusts` with respect to the ground truth clustering assignments in `truth`. 
# Return values
- `nbloss`: Normalised Binder loss (= 1 - rand index). Lower is better. 
- `ari`: Adjusted Rand index. Higher is better. 
- `vi`: Variation of Information distance (``\in [0, \log N]``). Lower is better.
- `nvi`: Normalised VI distance (``\in [0, 1]``).
- `id`: Information Distance (``\in [0, \log N]``). Lower is better. 
- `nid`: Normalised Information Distance (``\in [0, 1]``).
- `nmi`: Normalised Mutual Information (``\in [0, 1]``). Higher is better. 
"""
function evaluateclustering(clusts::ClustLabelVector, truth::ClustLabelVector)
    length(clusts) != length(truth) && throw(ArgumentError("Length of inputs must be equal."))
    n = length(clusts)
    ari, _, nbloss, _ = randindex(clusts, truth)
    vi = varinfo(clusts, truth)
    nvi = vi / log(n)
    id = infodist(clusts, truth; normalised = false)
    nid = id / log(n)
    nmi = mutualinfo(clusts, truth)
    return (nbloss = nbloss, ari = ari, vi = vi, nvi, id = id, nid = nid, nmi = nmi)
end

"""
    summarise(clusts::ClustLabelVector, 
    truth::ClustLabelVector, 
    printoutput = true) -> String
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
    if printoutput print(output) end
    return output
end