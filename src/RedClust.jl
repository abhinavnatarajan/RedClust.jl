@doc raw"""
# RedClust

## Introduction

RedClust is the main module of `RedClust.jl`, a Julia package for Bayesian clustering of high-dimensional Euclidean data using pairwise dissimilarities instead of the raw observations. 

Use `names(RedClust)` to get the export list of this module, and type `?name` to get help on a specific `name`. 

See https://abhinavnatarajan.github.io/RedClust.jl/ for detailed documentation.
"""
module RedClust

export 
# fit prior
fitprior, 
sampledist,
sampleK,

# MCMC sampler
runsampler,

# utility functions
adjacencymatrix, 
sortlabels, 
uppertriangle, 
generatemixture, 
makematrix, 

# summary functions
evaluateclustering, 
summarise,

# point estimate
getpointestimate,
binderloss,
infodist, 

# types
MCMCOptionsList, 
PriorHyperparamsList, 
MCMCData, 
MCMCResult,

# data
example_datasets, 
example_dataset

include("./types.jl")
include("./utils.jl")
include("./prior.jl")
include("./mcmc.jl")
include("./pointestimate.jl")
include("./summaries.jl")
include("./example_data.jl")

end