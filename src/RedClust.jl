@doc raw"""
# RedClust

## Introduction

[RedClust](https://github.com/abhinavnatarajan/RedClust.jl) is a [Julia](https://julialang.org/) package for Bayesian clustering of high-dimensional Euclidean data using pairwise dissimilarities instead of the raw observations. It uses an MCMC sampler to generate posterior samples from the space of all possible clustering structures on the data. 

## Installation
The package can be installed by typing `]add RedClust` into the Julia REPL or by the usual method:
```julia
using Pkg
Pkg.add("RedClust")
```

## Basic example
```julia
using RedClust
# Generate data
points, distM, clusts, probs, oracle_coclustering = 
	generatemixture(100, 10; α = 10, σ = 0.25, dim = 10)
# Let RedClust choose the best prior hyperparameters
params = fitprior(pnts, "k-means", false)
# Set the MCMC options
options = MCMCOptionsList(numiters = 5000)
data = MCMCData(points)
# Run the sampler
result = runsampler(data, options, params)
# Get a point estimate 
pointestimate, index = getpointestimate(result)
# Summary of point estimate
summarise(pointestimate, clusts)
```

## Citing this package
If you want to use this package in your work, please cite it as:

Natarajan, A., De Iorio, M., Heinecke, A., Mayer, E. and Glenn, S., 2022. ‘Cohesion and Repulsion in Bayesian Distance Clustering’, arXiv [2107.05414](https://arxiv.org/abs/2107.05414).
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
MCMCResult

# data
# example_datasets, 
# example_dataset

include("./types.jl")
include("./utils.jl")
include("./prior.jl")
include("./mcmc.jl")
include("./pointestimate.jl")
include("./summaries.jl")
include("./example_data.jl")

end