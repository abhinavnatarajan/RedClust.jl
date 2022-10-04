module RedClust

include("./Types.jl")
include("./Utility.jl")
include("./Prior.jl")
include("./MCMC.jl")
# include("./Plot.jl")
include("./SummaryFunctions.jl")


export fitprior, runsampler, # main functions
adjacencymatrix, sortlabels, uppertriangle, generatemixture, makematrix, pointestimate, # convenience
evaluateclustering, summarise, # summary functions 
MCMCOptionsList, PriorHyperparamsList, MCMCData, MCMCResult # types

end