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

include("./types.jl")
include("./utils.jl")
include("./prior.jl")
include("./mcmc.jl")
include("./pointestimate.jl")
# include("./Plot.jl")
include("./summaries.jl")


end