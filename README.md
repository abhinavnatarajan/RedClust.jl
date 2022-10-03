# RedClust

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://abhinavnatarajan.github.io/RedClust.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://abhinavnatarajan.github.io/RedClust.jl/dev/)
[![Build Status](https://github.com/abhinavnatarajan/RedClust.jl/actions/workflows/Test.yml/badge.svg?branch=master)](https://github.com/abhinavnatarajan/RedClust.jl/actions/workflows/Test.yml?query=branch%3Amaster)

## Introduction

[RedClust](https://github.com/abhinavnatarajan/RedClust.jl) is a [Julia](https://julialang.org/) package for Bayesian clustering of high-dimensional Euclidean data using pairwise dissimilarities instead of the raw observations. It uses an MCMC sampler to generate posterior samples from the space of all possible clustering structures on the data. 

## Installation
The package can be installed with 

	]add RedClust
or

	using Pkg
	Pkg.add("RedClust")
RedClust also requires [R](https://www.r-project.org/), and the R packages [salso](https://CRAN.R-project.org/package=salso) and [cluster](https://cran.r-project.org/package=cluster). If these are already installed, make sure the `R_HOME` environment variable is set to the `R` home directory. You could run `R.home()` in `R` to determine the location of this directory. If `R`, `salso`, and `cluster` are not found, they will be automatically installed by the Julia package [RCall.jl](https://github.com/JuliaInterop/RCall.jl). 

## Basic example

	using RedClust
	# Generate data
	pnts, distM, clusts, probs, oracle_coclustering = generatemixture(N, K; α = 10, σ = data_σ, dim = data_dim)
	# Let RedClust choose the best prior hyperparameters
	params = fitprior(pnts, "k-means", false).params
	# Set the MCMC options
	options = MCMCOptionsList(numiters = 5000)
	data = MCMCData(D = distM)
	# Run the sampler
	result = runsampler(data, options, params)

## Citing this package
If you want to use this package in your work, please cite it as:

Natarajan, A., De Iorio, M., Heinecke, A., Mayer, E. and Glenn, S., 2022. ‘Cohesion and Repulsion in Bayesian Distance Clustering’, arXiv [2107.05414](https://arxiv.org/abs/2107.05414).