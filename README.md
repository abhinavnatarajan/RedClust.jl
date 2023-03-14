# RedClust

[![GitHub Workflow Status (master)][github-CI-img]][github-CI-url]
[![License][license-img]][license-url]
[![Latest release][release-img]][release-url]
[![Code Coverage][codecov-img]][codecov-url]

## Documentation 

[![Development version documentation][docs-dev-img]][docs-dev-url]
[![Stable version documentation][docs-stable-img]][docs-stable-url]
[![arxiv paper link][arxiv-img]][arxiv-url]

Please see the detailed documentation above. 

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
A more elaborate example can be found in the detailed documentation. Examples from the paper and its supplementary material can be found in the ['examples' branch of this repository](https://github.com/abhinavnatarajan/RedClust.jl/tree/examples).

## Citing this package
If you want to use this package in your work, please cite it as:

Natarajan, A., De Iorio, M., Heinecke, A., Mayer, E. and Glenn, S. (2023). ‘Cohesion and Repulsion in Bayesian Distance Clustering’, _Journal of the Americal Statistical Association_. DOI: [10.1080/01621459.2023.2191821](https://doi.org/10.1080/01621459.2023.2191821).

For BibTeX users:

```latex
@article{NDI23,
  doi = {10.1080/01621459.2023.2191821},
  author = {Natarajan, Abhinav and De Iorio, Maria and Heinecke, Andreas and Mayer, Emanuel and Glenn, Simon},
  title = {Cohesion and Repulsion in Bayesian Distance Clustering},
  publisher = {Journal of the American Statistical Association},
  year = {2023}
}
```

[github-CI-img]: https://img.shields.io/github/actions/workflow/status/abhinavnatarajan/RedClust.jl/CI.yml?logo=github&labelColor=222222
[github-CI-url]: https://github.com/abhinavnatarajan/RedClust.jl/actions/workflows/CI.yml

[codecov-img]: https://img.shields.io/codecov/c/github/abhinavnatarajan/RedClust.jl?logo=codecov&labelColor=222222&logoColor=white
[codecov-url]: https://app.codecov.io/gh/abhinavnatarajan/RedClust.jl/

[release-img]: https://img.shields.io/github/v/release/abhinavnatarajan/RedClust.jl?display_name=tag&logo=SemVer&sort=semver&labelColor=222222
[release-url]: https://github.com/abhinavnatarajan/RedClust.jl/releases

[license-img]: https://img.shields.io/github/license/abhinavnatarajan/RedClust.jl?style=flat&labelColor=222222
[license-url]: https://github.com/abhinavnatarajan/RedClust.jl/blob/master/LICENSE

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue?style=flat&labelColor=222222
[docs-dev-url]: https://abhinavnatarajan.github.io/RedClust.jl/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue?style=flat&labelColor=222222
[docs-stable-url]: https://abhinavnatarajan.github.io/RedClust.jl/stable/

[arxiv-img]: https://img.shields.io/badge/arxiv-2107.05414-red?style=flat&labelColor=222222
[arxiv-url]: https://arxiv.org/abs/2107.05414
