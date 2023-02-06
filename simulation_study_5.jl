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

# Folder to save things to
saveloc = joinpath(@__DIR__, "figures", "simulation_study_5")
if !isdir(saveloc)
    mkpath(saveloc)
end

# Get the data
begin
    data = example_dataset(dataset_number)
    points, distmatrix, clusts, probs, oracle_coclustering = data
    N = length(points)
end

# RedClust includes the function `fitprior2` to heuristically choose prior hyperparameters based on the data. This is a variation on `fitprior`, using an ensemble of possible clusterings to estimate the prior density on distances.
params = fitprior2(points, "k-means", false)

# We can check how good the chosen prior hyperparameters are by comparing the empirical distribution of distances to the (marginal) prior predictive distribution.
begin
    empirical_intracluster = uppertriangle(distmatrix)[
        uppertriangle(adjacencymatrix(clusts)) .== 1]
    empirical_intercluster = uppertriangle(distmatrix)[
        uppertriangle(adjacencymatrix(clusts)) .== 0]
        pred_intracluster = sampledist(params, "intracluster", 10000)
        pred_intercluster = sampledist(params, "intercluster", 10000)
        density(pred_intercluster,
        label="Simulated ICD",
        linewidth = 2, linestyle = :dash, legend=false)
        density!(empirical_intercluster,
        label="Empirical ICD",
        linewidth = 2, color = 1)
        density!(pred_intracluster,
        label="Simulated WCD", xlabel = "Distance", ylabel = "Density",
        linewidth = 2, linestyle = :dash, color = 2)
        density!(empirical_intracluster,
        label="Empirical WCD",
        color = 2, linewidth = 2)
        title!("Distances: Prior Predictive vs Empirical Distribution")
end
savefig(joinpath(saveloc,"predicted_density_distances.pdf"))

# We can also evaluate the prior hyperparameters by checking the marginal predictive distribution on ``K`` (the number of clusters).
begin
    Ksamples = sampleK(params, 10000, N)
    histogram_pmf(Ksamples, binwidth = 3, legend = false,
    xticks=collect(0:20:maximum(Ksamples)),
    yticks=collect(0:0.04:1.0),
    xlabel = "\$K\$", ylabel = "Probability", 
    title = "Marginal Prior Predictive Distribution of \$K\$",
    label="Empirical density")
end
savefig(joinpath(saveloc,"prior_predictive_K.pdf"))

# Running the MCMC is straightforward. We set up the MCMC options using `MCMCOptionsList`.
options = MCMCOptionsList(numiters=numiters)

# We then set up the input data using `MCMCData`.
data = MCMCData(points)

# We can then run the sampler using `runsampler`.
result = runsampler(data, options, params)

# Plot the posterior coclustering matrix
sqmatrixplot(result.posterior_coclustering,
title="Posterior Coclustering Probabilities")
savefig(joinpath(saveloc,"posterior_coclustering_matrix.pdf"))

# Plot the posterior distribution of K:
histogram_pmf(result.K,
xlabel = "\$K\$", ylabel = "Probability", title = "Posterior Distribution of \$K\$")
savefig(joinpath(saveloc,"K_posterior.pdf"))

# Plot the posterior distribution of r:
begin
    histogram(result.r, normalize = :pdf,
    legend=true, 
    label="Empirical density", ylabel = "Density", xlabel = "\$r\$", 
    title = "Posterior Distribution of \$r\$")
    density!(result.r, 
    color=:black, linewidth = 2, linestyle=:dash, 
    label="Kernel estimate")
end
savefig(joinpath(saveloc,"r_posterior.pdf"))

# Plot the posterior distribution of p:
begin
    histogram(result.p, normalize = :pdf,
    ylabel = "Density", xlabel = "\$p\$",
    title = "Posterior Distribution of \$p\$",
    label = "Empirical density",
    legend=true)
    density!(result.p, color=:black, linewidth = 2, linestyle=:dash,
    label = "Kernel estimate", legend_position = :topleft)
end
savefig(joinpath(saveloc,"p_posterior.pdf"))

# Plot the traceplot of the autocorrelation function of ``K``:
plot(result.K_acf, legend = false, linewidth = 1,
xlabel = "Lag", ylabel = "Autocorrelation",
title = "Autocorrelation Function of \$K\$")
savefig(joinpath(saveloc,"K_acf.pdf"))

# Plot the traceplot of the autocorrelation function of ``r``:
plot(result.r_acf, legend = false, linewidth = 1,
xlabel = "Lag", ylabel = "Autocorrelation",
title = "Autocorrelation Function of \$r\$")
savefig(joinpath(saveloc,"r_acf.pdf"))

# Plot the traceplot of the autocorrelation function of ``p``:
plot(result.p_acf, legend = false, linewidth = 1,
xlabel = "Lag", ylabel = "Autocorrelation",
title = "Autocorrelation Function of \$p\$")
savefig(joinpath(saveloc,"p_acf.pdf"))

# Check the trace plot of the log-likelihood to make sure the MCMC is moving well:
plot(result.loglik, legend = false, linewidth = 1,
xlabel = "Iteration", ylabel = "Log likelihood",
xformatter=(x -> (x == 0) ? "0" : string(Int(floor(x/10000))) * "e4"),
title = "Log-Likelihood Trace Plot")
savefig(joinpath(saveloc,"loglik.pdf"))

# Check the trace plot of the log-posterior:
plot(result.logposterior, legend = false, linewidth = 1,
xlabel = "Iteration", ylabel = "Log posterior",
title = "Log-Posterior Trace Plot",
xformatter=(x -> (x == 0) ? "0" : string(Int(floor(x/10000))) * "e4"))
savefig(joinpath(saveloc,"logposterior.pdf"))

# Use the salso algorithm to find a point estimate
result_z = permutedims(reduce(hcat, result.clusts))
pointestimate = rcopy(R"""
    library(salso)
    salso($result_z)
    """)
# Visualise the adjacency matrix of the point-estimate
sqmatrixplot(adjacencymatrix(pointestimate), 
title = "Adjacency Matrix of Point Estimate", colorbar=false)
savefig(joinpath(saveloc,"predicted_coclustering_matrix.pdf"))

# We can check the accuracy of the point estimate in terms of clustering metrics.
summarise(pointestimate, clusts)
# Save the results into a file
begin
local summary_file = open(joinpath(saveloc, "summary.txt"), "w")
show(summary_file, "text/plain", result)
println(summary_file, "")
summarise(summary_file, pointestimate, clusts)
close(summary_file)
end