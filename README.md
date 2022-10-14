# RedClust

![GitHub Workflow Status (master)](https://img.shields.io/github/workflow/status/abhinavnatarajan/RedClust.jl/Run%20tests/master?style=flat&logo=github)
[![License](https://img.shields.io/github/license/abhinavnatarajan/RedClust.jl?style=flat)](https://github.com/abhinavnatarajan/RedClust.jl/blob/master/LICENSE)
[![GitHub tag (latest SemVer)](https://img.shields.io/github/v/tag/abhinavnatarajan/RedClust.jl?style=flat) ](https://github.com/abhinavnatarajan/RedClust.jl/releases)

![GitHub Documentation Status (latest)](https://img.shields.io/github/workflow/status/abhinavnatarajan/RedClust.jl/Documentation?label=Documentation)
[![Stable](https://img.shields.io/badge/docs-stable-blue?style=flat)](https://abhinavnatarajan.github.io/RedClust.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue?style=flat)](https://abhinavnatarajan.github.io/RedClust.jl/dev/)
[![arXiv](https://img.shields.io/badge/arxiv-2107.05414-red)](https://arxiv.org/abs/2107.05414)

## Introduction

[RedClust](https://github.com/abhinavnatarajan/RedClust.jl) is a [Julia](https://julialang.org/) package for Bayesian clustering of high-dimensional Euclidean data using pairwise dissimilarities instead of the raw observations. It uses an MCMC sampler to generate posterior samples from the space of all possible clustering structures on the data. 

## Installation
The package can be installed by typing `]add RedClust` into the Julia REPL or by the usual method:
```julia
using Pkg
Pkg.add("RedClust")
```
RedClust also requires [R](https://www.r-project.org/) and the R package [`salso`](https://CRAN.R-project.org/package=salso). If R is already installed, make sure the `R_HOME` environment variable is set to the R home directory (you could run `R.home()` in R to determine the location of this directory). If R or `salso` are not found, they are automatically installed during package installation.   

## Basic example
```julia
using RedClust
# Generate data
points, distM, clusts, probs, oracle_coclustering = 
	generatemixture(N, K; α = 10, σ = data_σ, dim = data_dim)
# Let RedClust choose the best prior hyperparameters
params = fitprior(pnts, "k-means", false)
# Set the MCMC options
options = MCMCOptionsList(numiters = 5000)
data = MCMCData(points)
# Run the sampler
result = runsampler(data, options, params)
# Get a point estimate 
pointestimate, index = getpointestimate(result)
# Summary of MCMC and point estimate
summarise(result)
summarise(pointestimate, clusts)
```

## Citing this package
If you want to use this package in your work, please cite it as:

Natarajan, A., De Iorio, M., Heinecke, A., Mayer, E. and Glenn, S., 2022. ‘Cohesion and Repulsion in Bayesian Distance Clustering’, arXiv [2107.05414](https://arxiv.org/abs/2107.05414).
