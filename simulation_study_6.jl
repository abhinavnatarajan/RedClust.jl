# We begin by setting up the necessary includes.
using Pkg
Pkg.activate(@__DIR__)
# Pkg.develop(Pkg.PackageSpec(path = joinpath(@__DIR__, "..")))
Pkg.instantiate()
using RCall: rcopy, @R_str
using RedClust, Plots, StatsPlots, BayesianMixtures
using Random: seed!
using StatsBase: counts, mean
using LinearAlgebra: triu, diagind
using Clustering: kmeans, kmedoids
using FileIO: open, close
include("utils_for_examples.jl")

dataset_number = 1
numiters = 50000 # global for all algorithms

# seed the default RNG for stability
seed!(44)

# Get the data
begin
    data = example_dataset(dataset_number)
    points, distmatrix, clusts, probs, oracle_coclustering = data
    N = length(points)
end

# RedClust includes the function `fitprior` to heuristically choose prior hyperparameters based on the data.
params = fitprior(points, "k-means", false)
params.repulsion = false

# Running the MCMC is straightforward. We set up the MCMC options using `MCMCOptionsList`.
options = MCMCOptionsList(numiters=numiters)

# We then set up the input data using `MCMCData`.
data = MCMCData(points)

result = MCMCResult[]

for (idx, maxK) in enumerate([10, 25, 50, 75, 100])
    # Folder to save things to
    saveloc = joinpath(@__DIR__, "figures", "simulation_study_6", "maxK" * string(maxK))
    if !isdir(saveloc)
        mkpath(saveloc)
    end
    params.maxK = maxK
    # We can then run the sampler using `runsampler`.
    # push!(result, runsampler(data, options, params))

    # Plot the posterior coclustering matrix
    sqmatrixplot(result[idx].posterior_coclustering,
    title="Posterior Coclustering Probabilities")
    savefig(joinpath(saveloc, "posterior_coclustering_matrix.pdf"))

    # Plot the posterior distribution of K:
    temp = [1, 1, 2, 3, 3]
    histogram_pmf(result[idx].K, 
    xticks=minimum(result[idx].K):temp[idx]:maximum(result[idx].K),
    xlabel = "\$K\$", ylabel = "Probability", 
    title = "Posterior Distribution of \$K\$ (No repulsion)")
    savefig(joinpath(saveloc, "K_posterior.pdf"))

    # Plot the posterior distribution of r:
    begin
        histogram(result[idx].r, normalize = :pdf,
        legend=true, 
        label="Empirical density", ylabel = "Density", xlabel = "\$r\$", 
        title = "Posterior Distribution of \$r\$")
        density!(result[idx].r, 
        color=:black, linewidth = 2, linestyle=:dash, 
        label="Kernel estimate")
    end
    savefig(joinpath(saveloc,"r_posterior.pdf"))

    # Plot the posterior distribution of p:
    begin
        histogram(result[idx].p, normalize = :pdf,
        ylabel = "Density", xlabel = "\$p\$",
        title = "Posterior Distribution of \$p\$",
        label = "Empirical density",
        legend=false)
        density!(result[idx].p, color=:black, linewidth = 2, linestyle=:dash,
        label = "Kernel estimate")
    end
    savefig(joinpath(saveloc,"p_posterior.pdf"))
end