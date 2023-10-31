@doc raw"""
# RedClust

## Introduction

RedClust is the main module of `RedClust.jl`, a Julia package for Bayesian clustering of high-dimensional Euclidean data using pairwise dissimilarities instead of the raw observations.

Use `names(RedClust)` to get the export list of this module, and type `?name` to get help on a specific `name`.

See https://abhinavnatarajan.github.io/RedClust.jl/ for detailed documentation.
"""
module RedClust

# mcmc.jl
using Clustering: kmeans, kmedoids, mutualinfo, randindex, varinfo, vmeasure
using Dates
using Distances: Euclidean, pairwise
using Distributions: Beta, Dirichlet, Distributions, Gamma, MvNormal, Normal, Uniform,
    fit_mle, logpdf, pdf, rate, shape, truncated
using HDF5: close, create_group, h5open
using LinearAlgebra: Diagonal, I
using LoopVectorization: @turbo
using Printf: @sprintf
using ProgressBars: ProgressBar
using Random: AbstractRNG, TaskLocalRNG, rand, seed!
using SpecialFunctions: logbeta, loggamma
using StaticArrays
using StatsBase:
    autocor, counts, entropy, levelsmap, mean, mean_and_var, sample, std, var, wsample


export
    # fit prior
    fitprior,
    fitprior2,
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
