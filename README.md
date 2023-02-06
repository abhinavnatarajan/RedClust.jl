## RedClust.jl Examples
This branch of the repository contains the code used to run the simulations and generate the figures displayed in the following paper:

[![arxiv paper link][arxiv-img]][arxiv-url] Natarajan, A., De Iorio, M., Heinecke, A., Mayer, E. and Glenn, S., 2021. ‘Cohesion and Repulsion in Bayesian Distance Clustering’, arXiv [2107.05414](https://arxiv.org/abs/2107.05414).


## Instructions for Reproducing the Results in the Paper
In order to reproduce the figures from the paper, you will need an environment with the following:
- Julia 1.8.2
- R 4.2.2
- The R package `salso`, available on CRAN.
Before running the examples, please modify and run the file `setup.jl` in order to ensure that dependencies are correctly set up. It is possible that your results may differ slightly from the ones in the paper due to variability in the random number generator on your system, seeding, and versions of the various Julia packages and their dependencies. 

[arxiv-img]: https://img.shields.io/badge/arxiv-2107.05414-red?style=flat&labelColor=222222
[arxiv-url]: https://arxiv.org/abs/2107.05414

