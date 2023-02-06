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
saveloc = joinpath(@__DIR__, "figures", "simulation_study_4")
if !isdir(saveloc)
    mkpath(saveloc)
end

# Get the data
begin
    data = example_dataset(dataset_number)
    points, distmatrix, clusts, probs, oracle_coclustering = data
    N = length(points)
end

# RedClust includes the function `fitprior` to heuristically choose prior hyperparameters based on the data. Here we test the robustness of the prior parameters estimated by varying the value of ``K`` chosen by the elbow method.
params = fitprior(points, "k-means", false)
K_elbow = params.K_initial
params_Kp2 = fitprior(points, "k-means", false; Kmin = K_elbow+3, Kmax = K_elbow+3)
params_Km2 = fitprior(points, "k-means", false; Kmin = K_elbow-3, Kmax = K_elbow-3)

# We can check how good the chosen prior hyperparameters are by comparing the empirical distribution of distances to the (marginal) prior predictive distribution.
begin
    empirical_intracluster = uppertriangle(distmatrix)[
        uppertriangle(adjacencymatrix(clusts)) .== 1]
    empirical_intercluster = uppertriangle(distmatrix)[
        uppertriangle(adjacencymatrix(clusts)) .== 0]
    pred_intracluster = sampledist(params_Kp2, "intracluster", 10000)
    pred_intercluster = sampledist(params_Kp2, "intercluster", 10000)
    density(pred_intercluster,
    label="Simulated ICD", xlabel = "Distance", ylabel = "Density",
    linewidth = 2, linestyle = :dash, legend=false)
    density!(empirical_intercluster,
    label="Empirical ICD",
    color = 1, linewidth = 2)
    density!(pred_intracluster,
    label="Simulated WCD",
    linewidth = 2, linestyle = :dash, color = 2)
    density!(empirical_intracluster,
    label="Empirical WCD",
    linewidth = 2, color = 2)
    title!("Distances: Prior Predictive vs Empirical Distribution")
end
savefig(joinpath(saveloc,"predicted_density_distances_Kp2.pdf"))
begin
    pred_intracluster = sampledist(params_Km2, "intracluster", 10000)
    pred_intercluster = sampledist(params_Km2, "intercluster", 10000)
    density(pred_intercluster,
    label="Simulated ICD", xlabel = "Distance", ylabel = "Density",
    linewidth = 2, linestyle = :dash, legend=false)
    density!(empirical_intercluster,
    label="Empirical ICD",
    color = 1, linewidth = 2)
    density!(pred_intracluster,
    label="Simulated WCD",
    linewidth = 2, linestyle = :dash, color = 2)
    density!(empirical_intracluster,
    label="Empirical WCD",
    linewidth = 2, color = 2)
    title!("Distances: Prior Predictive vs Empirical Distribution")
end
savefig(joinpath(saveloc,"predicted_density_distances_Km2.pdf"))

# We can also evaluate the prior hyperparameters by checking the marginal predictive distribution on ``K`` (the number of clusters).
begin
    Ksamples_Kp2 = sampleK(params_Kp2, 10000, N)
    histogram_pmf(Ksamples_Kp2, binwidth = 3, legend = false,
    xticks=collect(0:20:maximum(Ksamples_Kp2)),
    yticks=collect(0:0.04:1.0),
    xlabel = "\$K\$", ylabel = "Probability", 
    title = "Marginal Prior Predictive Distribution of \$K\$",
    label="Empirical density")
end
savefig(joinpath(saveloc,"prior_predictive_K_Kp2.pdf"))
begin
    Ksamples_Km2 = sampleK(params_Km2, 10000, N)
    histogram_pmf(Ksamples_Km2, binwidth = 3, legend = false,
    xticks=collect(0:20:maximum(Ksamples_Km2)),
    yticks=collect(0:0.04:1.0),
    xlabel = "\$K\$", ylabel = "Probability", 
    title = "Marginal Prior Predictive Distribution of \$K\$",
    label="Empirical density")
end
savefig(joinpath(saveloc,"prior_predictive_K_Km2.pdf"))

# Running the MCMC is straightforward. We set up the MCMC options using `MCMCOptionsList`.
options = MCMCOptionsList(numiters=numiters)

# We then set up the input data using `MCMCData`.
data = MCMCData(points)

# We can then run the sampler using `runsampler`.
result_Kp2 = runsampler(data, options, params_Kp2)
result_Km2 = runsampler(data, options, params_Km2)

# Plot the posterior coclustering matrix
sqmatrixplot(result_Kp2.posterior_coclustering,
title="Posterior Coclustering Probabilities")
savefig(joinpath(saveloc,"posterior_coclustering_matrix_Kp2.pdf"))
# Plot the posterior coclustering matrix
sqmatrixplot(result_Km2.posterior_coclustering,
title="Posterior Coclustering Probabilities")
savefig(joinpath(saveloc,"posterior_coclustering_matrix_Km2.pdf"))

# Plot the posterior distribution of K:
histogram_pmf(result_Kp2.K,
xlabel = "\$K\$", ylabel = "Probability", title = "Posterior Distribution of \$K\$")
savefig(joinpath(saveloc,"K_posterior_Kp2.pdf"))
histogram_pmf(result_Km2.K,
xlabel = "\$K\$", ylabel = "Probability", title = "Posterior Distribution of \$K\$")
savefig(joinpath(saveloc,"K_posterior_Km2.pdf"))

# Plot the posterior distribution of r:
begin
    histogram(result_Kp2.r, normalize = :pdf,
    legend=true, 
    label="Empirical density", ylabel = "Density", xlabel = "\$r\$", 
    title = "Posterior Distribution of \$r\$")
    density!(result_Kp2.r, 
    color=:black, linewidth = 2, linestyle=:dash, 
    label="Kernel estimate")
end
savefig(joinpath(saveloc,"r_posterior_Kp2.pdf"))
begin
    histogram(result_Km2.r, normalize = :pdf,
    legend=true, 
    label="Empirical density", ylabel = "Density", xlabel = "\$r\$", 
    title = "Posterior Distribution of \$r\$")
    density!(result_Km2.r, 
    color=:black, linewidth = 2, linestyle=:dash, 
    label="Kernel estimate")
end
savefig(joinpath(saveloc,"r_posterior_Km2.pdf"))

# Plot the posterior distribution of p:
begin
    histogram(result_Kp2.p, normalize = :pdf,
    ylabel = "Density", xlabel = "\$p\$",
    title = "Posterior Distribution of \$p\$",
    label = "Empirical density")
    density!(result_Kp2.p, color=:black, linewidth = 2, linestyle=:dash,
    label = "Kernel estimate",
    legend_position = :topleft)
end
savefig(joinpath(saveloc,"p_posterior_Kp2.pdf"))
begin
    histogram(result_Km2.p, normalize = :pdf,
    ylabel = "Density", xlabel = "\$p\$",
    title = "Posterior Distribution of \$p\$",
    label = "Empirical density")
    density!(result_Km2.p, 
    color=:black, linewidth = 2, linestyle=:dash,
    label = "Kernel estimate",
    legend_position = :topleft)
end
savefig(joinpath(saveloc,"p_posterior_Km2.pdf"))

# Plot the traceplot of the autocorrelation function of ``K``:
plot(result_Kp2.K_acf, legend = false, linewidth = 1,
xlabel = "Lag", ylabel = "Autocorrelation",
title = "Autocorrelation Function of \$K\$")
savefig(joinpath(saveloc,"K_acf_Kp2.pdf"))
plot(result_Km2.K_acf, legend = false, linewidth = 1,
xlabel = "Lag", ylabel = "Autocorrelation",
title = "Autocorrelation Function of \$K\$")
savefig(joinpath(saveloc,"K_acf_Km2.pdf"))

# Plot the traceplot of the autocorrelation function of ``r``:
plot(result_Kp2.r_acf, legend = false, linewidth = 1,
xlabel = "Lag", ylabel = "Autocorrelation",
title = "Autocorrelation Function of \$r\$")
savefig(joinpath(saveloc,"r_acf_Kp2.pdf"))
plot(result_Km2.r_acf, legend = false, linewidth = 1,
xlabel = "Lag", ylabel = "Autocorrelation",
title = "Autocorrelation Function of \$r\$")
savefig(joinpath(saveloc,"r_acf_Km2.pdf"))

# Plot the traceplot of the autocorrelation function of ``p``:
plot(result_Kp2.p_acf, legend = false, linewidth = 1,
xlabel = "Lag", ylabel = "Autocorrelation",
title = "Autocorrelation Function of \$p\$")
savefig(joinpath(saveloc,"p_acf_Kp2.pdf"))
plot(result_Km2.p_acf, legend = false, linewidth = 1,
xlabel = "Lag", ylabel = "Autocorrelation",
title = "Autocorrelation Function of \$p\$")
savefig(joinpath(saveloc,"p_acf_Km2.pdf"))

# Check the trace plot of the log-likelihood to make sure the MCMC is moving well:
plot(result_Kp2.loglik, legend = false, linewidth = 1,
xlabel = "Iteration", ylabel = "Log likelihood",
xformatter=(x -> (x == 0) ? "0" : string(Int(floor(x/10000))) * "e4"),
title = "Log-Likelihood Trace Plot")
savefig(joinpath(saveloc,"loglik_Kp2.pdf"))
plot(result_Km2.loglik, legend = false, linewidth = 1,
xlabel = "Iteration", ylabel = "Log likelihood",
xformatter=(x -> (x == 0) ? "0" : string(Int(floor(x/10000))) * "e4"),
title = "Log-Likelihood Trace Plot")
savefig(joinpath(saveloc,"loglik_Km2.pdf"))

# Check the trace plot of the log-posterior:
plot(result_Kp2.logposterior, legend = false, linewidth = 1,
xlabel = "Iteration", ylabel = "Log posterior",
xformatter=(x -> (x == 0) ? "0" : string(Int(floor(x/10000))) * "e4"),
title = "Log-Posterior Trace Plot")
savefig(joinpath(saveloc,"logposterior_Kp2.pdf"))
plot(result_Km2.logposterior, legend = false, linewidth = 1,
xlabel = "Iteration", ylabel = "Log posterior",
xformatter=(x -> (x == 0) ? "0" : string(Int(floor(x/10000))) * "e4"),
title = "Log-Posterior Trace Plot")
savefig(joinpath(saveloc,"logposterior_Km2.pdf"))

# Use the salso algorithm to find a point estimate
result_z_Kp2 = permutedims(reduce(hcat, result_Kp2.clusts))
pointestimate_Kp2 = rcopy(R"""
    library(salso)
    salso($result_z_Kp2)
    """)
# Visualise the adjacency matrix of the point-estimate
sqmatrixplot(adjacencymatrix(pointestimate_Kp2), 
title = "Adjacency Matrix of Point Estimate", colorbar=false)
savefig(joinpath(saveloc,"predicted_coclustering_matrix_Kp2.pdf"))
result_z_Km2 = permutedims(reduce(hcat, result_Km2.clusts))
pointestimate_Km2 = rcopy(R"""
    library(salso)
    salso($result_z_Km2)
    """)
# Visualise the adjacency matrix of the point-estimate
sqmatrixplot(adjacencymatrix(pointestimate_Km2), 
title = "Adjacency Matrix of Point Estimate", colorbar=false)
savefig(joinpath(saveloc,"predicted_coclustering_matrix_Km2.pdf"))


# We can check the accuracy of the point estimate in terms of clustering metrics.
summarise(pointestimate_Kp2, clusts)
summarise(pointestimate_Km2, clusts)
# Save the results into a file
begin
local summary_file = open(joinpath(saveloc, "summary_Kp2.txt"), "w")
show(summary_file, "text/plain", result_Kp2)
println(summary_file, "")
summarise(summary_file, pointestimate_Kp2, clusts)
close(summary_file)
end
begin
    local summary_file = open(joinpath(saveloc, "summary_Km2.txt"), "w")
    show(summary_file, "text/plain", result_Km2)
    println(summary_file, "")
    summarise(summary_file, pointestimate_Km2, clusts)
    close(summary_file)
end