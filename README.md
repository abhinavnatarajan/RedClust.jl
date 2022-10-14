# RedClust
*A Julia package for Bayesian clustering of high-dimensional data using pairwise dissimilarities and repulsion.*

[![GitHub Workflow Status (master)][github-CI-img]][github-CI-url]
[![License][license-img]][license-url]
[![Latest release][release-img]][release-url]
[![Code Coverage][codecov-img]][codecov-url]

[![Development version documentation][docs-stable-img]][docs-stable-url]
[![Stable version documentation][docs-stable-img]][docs-stable-url]
[![arxiv paper link][arxiv-img]][arxiv-url]

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

[github-CI-img]: https://img.shields.io/github/workflow/status/abhinavnatarajan/RedClust.jl/CI?label=CI&logo=github
[github-CI-url]: https://github.com/abhinavnatarajan/RedClust.jl/actions/workflows/CI.yml

[codecov-img]: https://img.shields.io/codecov/c/github/abhinavnatarajan/RedClust.jl?logo=codecov
[codecov-url]: https://app.codecov.io/gh/abhinavnatarajan/RedClust.jl/

[release-img]: https://img.shields.io/github/v/release/abhinavnatarajan/RedClust.jl?display_name=tag&logo=SemVer&sort=semver
[release-url]: https://github.com/abhinavnatarajan/RedClust.jl/releases

[license-img]: https://img.shields.io/github/license/abhinavnatarajan/RedClust.jl?style=flat
[license-url]: https://github.com/abhinavnatarajan/RedClust.jl/blob/master/LICENSE

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue?style=flat
[docs-dev-url]: https://abhinavnatarajan.github.io/RedClust.jl/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue?style=flat
[docs-stable-url]: https://abhinavnatarajan.github.io/RedClust.jl/stable/

[arxiv-img]: https://img.shields.io/badge/arxiv-2107.05414-red
[arxiv-url]: https://arxiv.org/abs/2107.05414