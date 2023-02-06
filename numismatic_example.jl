# We begin by setting up the necessary includes.
using Pkg
Pkg.activate(@__DIR__)
# Pkg.develop(Pkg.PackageSpec(path = joinpath(@__DIR__, "..")))
Pkg.instantiate()
using RCall: rcopy, @R_str
using RedClust, Plots, StatsPlots, BayesianMixtures 
import CSV: File
using DataFrames: DataFrame, nrow, ncol
using Random: seed!
using StatsBase: counts, mean
using LinearAlgebra: triu, diagind
using Clustering: kmeans, kmedoids
using Distances: Euclidean, pairwise
using MultivariateStats: fit, MetricMDS, predict
using FileIO: open, close
include("utils_for_examples.jl")

numiters = 50000 # global for all algorithms

seed!(44)

# Folder to save things to
saveloc = joinpath(@__DIR__, "figures", "numismatic_example")
if !isdir(saveloc)
    mkpath(saveloc)
end
dataloc = joinpath(@__DIR__, "coins_data")

# Get the raw and processed distances
begin
    local rawdistancedf = DataFrame(File(joinpath(dataloc, "coins_raw_distances.csv"), header = true))
    local temp = Bool.([rawdistancedf[i, 2] for i in 1:nrow(rawdistancedf)])
    within_cluster_raw_distances = rawdistancedf[temp, 1]
    inter_cluster_raw_distances = rawdistancedf[.!temp, 1]
    histogram(inter_cluster_raw_distances, label = "ICD", 
    title = "Histogram of Raw Distances", 
    ylabel = "Frequency", xlabel = "Distance", 
    legend=false, size=(600/400*470/2.5, 470/2.5))
    histogram!(within_cluster_raw_distances, label = "WCD")
end
savefig(joinpath(saveloc, "histogram_raw_distances.pdf"))
begin
    local distancedf = DataFrame(File(joinpath(dataloc, "coins_distance_matrix.csv"), header = false))
    N = nrow(distancedf)
    distmatrix = [distancedf[i, j] for i in 1:nrow(distancedf), j in 1:ncol(distancedf)]
    local maxfinitedist = maximum(distmatrix[.!(distmatrix .== Inf)])
    local minfinitedist = minimum(distmatrix[distmatrix .> 0])
    distmatrix[distmatrix .== Inf] .= maxfinitedist * 1.1
    distmatrix[distmatrix .== -Inf] .= minfinitedist * 0.9
    distmatrix[[CartesianIndex(i, i) for i in 1:N]] .= 0
    clusts = DataFrame(File(joinpath(dataloc, "coins_labels.csv"))).Group
    K = length(unique(clusts))
    local mdsmodel = fit(MetricMDS, distmatrix;
    distances = true, maxiter = 1000)
    points = predict(mdsmodel)
    MDSdistmatrix = pairwise(Euclidean(), points, dims = 2)
    points = [points[:, i] for i in axes(points, 2)]
end

# We can visualise the true adjacency matrix of the observations with respect to the true clusters that they were drawn from. 
sqmatrixplot(adjacencymatrix(clusts), title = "True adjacency matrix", colorbar=false)
savefig(joinpath(saveloc,"true_adjacency_matrix.pdf"))

# We can visualise the matrix of pairwise distances between the observations.
sqmatrixplot(distmatrix, title = "Matrix of Pairwise Distances")
savefig(joinpath(saveloc,"pairwise_distances_matrix.pdf"))

# We can also plot the histogram of distances, grouped by whether they are inter-cluster distances (ICD) or within-cluster distances (WCD).
begin
    empirical_intracluster = uppertriangle(distmatrix)[
        uppertriangle(adjacencymatrix(clusts)) .== 1]
    empirical_intercluster = uppertriangle(distmatrix)[
        uppertriangle(adjacencymatrix(clusts)) .== 0]
    histogram(empirical_intercluster,
    label="ICD", xlabel = "Distance", ylabel="Frequency",
    title = "Observed distribution of distances", 
    xticks = collect(0:0.015:0.06),
    bins = minimum(empirical_intercluster):0.002:maximum(empirical_intercluster),
    legend=false, size=(600/400*470/2.5, 470/2.5))
    histogram!(empirical_intracluster,
    bins = minimum(empirical_intracluster):0.002:maximum(empirical_intracluster),
    label="WCD")
end
savefig(joinpath(saveloc,"histogram_distances.pdf"))

# We can also plot the histogram of distances after MDS embedding
begin
    empirical_intracluster = uppertriangle(MDSdistmatrix)[
        uppertriangle(adjacencymatrix(clusts)) .== 1]
    empirical_intercluster = uppertriangle(MDSdistmatrix)[
        uppertriangle(adjacencymatrix(clusts)) .== 0]
    histogram(empirical_intercluster,
    label="ICD", xlabel = "Distance", ylabel="Frequency",
    title = "Distribution of distances after MDS", 
    bins = minimum(empirical_intercluster):0.002:maximum(empirical_intercluster),
    legend=false, size=(600/400*470/2.5, 470/2.5))
    histogram!(empirical_intracluster,
    bins = minimum(empirical_intracluster):0.002:maximum(empirical_intracluster),
    label="WCD")
end
savefig(joinpath(saveloc,"histogram_distances_MDS.pdf"))

# RedClust includes the function `fitprior` to heuristically choose prior hyperparameters based on the data.
params = fitprior(distmatrix, "k-medoids", true)

# We can check how good the chosen prior hyperparameters are by comparing the empirical distribution of distances to the (marginal) prior predictive distribution.
begin
    empirical_intracluster = uppertriangle(distmatrix)[
        uppertriangle(adjacencymatrix(clusts)) .== 1]
    empirical_intercluster = uppertriangle(distmatrix)[
        uppertriangle(adjacencymatrix(clusts)) .== 0]
    pred_intracluster = sampledist(params, "intracluster", 10000)
    pred_intercluster = sampledist(params, "intercluster", 10000)
    density(pred_intercluster,
    label="Simulated ICD", xlabel = "Distance", ylabel = "Density",
    linewidth = 2, linestyle = :dash, legend = false)
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
savefig(joinpath(saveloc,"predicted_density_distances.pdf"))

# We can also evaluate the prior hyperparameters by checking the marginal predictive distribution on ``K`` (the number of clusters).
begin
    Ksamples = sampleK(params, 10000, N)
    histogram_pmf(Ksamples, legend = false,
    xticks=collect(0:10:maximum(Ksamples)),
    xlabel = "\$K\$", ylabel = "Probability", title = "Marginal Prior Predictive Distribution of \$K\$")
end
savefig(joinpath(saveloc,"prior_predictive_K.pdf"))

# Running the MCMC is straightforward. We set up the MCMC options using `MCMCOptionsList`.
options = MCMCOptionsList(numiters=numiters)

# We then set up the input data using `MCMCData`.
data = MCMCData(distmatrix)

# We can then run the sampler using `runsampler`.
result = runsampler(data, options, params)

# Plot the posterior coclustering matrix
sqmatrixplot(result.posterior_coclustering,
title="Posterior Coclustering Probabilities")
savefig(joinpath(saveloc,"posterior_coclustering_matrix.pdf"))

# Plot the posterior distribution of K:
histogram_pmf(result.K,
xlabel = "\$K\$", ylabel = "Probability", 
title = "Posterior Distribution of \$K\$",
xticks = minimum(result.K):2:maximum(result.K))
savefig(joinpath(saveloc,"K_posterior.pdf"))

# Plot the posterior distribution of r:
begin
    histogram(result.r, normalize = :pdf, 
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
    legend_position = :topleft)
    density!(result.p, color=:black, linewidth = 2, 
    linestyle=:dash,
    label = "Kernel estimate")
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
title = "Log-Likelihood Trace Plot")
savefig(joinpath(saveloc,"loglik.pdf"))

# Check the trace plot of the log-posterior:
plot(result.logposterior, legend = false, linewidth = 1,
xlabel = "Iteration", ylabel = "Log posterior",
title = "Log-Posterior Trace Plot")
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

# Next we run the model without repulsion
params_norepulsion = deepcopy(params)
params_norepulsion.repulsion = false
result_norepulsion = runsampler(data, options, params_norepulsion)

# Plot the posterior coclustering matrix
sqmatrixplot(result_norepulsion.posterior_coclustering,
title="Posterior Coclustering Probabilities (No Repulsion)")
savefig(joinpath(saveloc,"posterior_coclustering_matrix_norepulsion.pdf"))

# Plot the posterior distribution of K:
histogram_pmf(result_norepulsion.K,
xlabel = "\$K\$", ylabel = "Probability", title = "Posterior Distribution of \$K\$ (No repulsion)")
savefig(joinpath(saveloc,"K_posterior_norepulsion.pdf"))

# Use the salso algorithm to find a point estimate
# result_z_norepulsion = permutedims(reduce(hcat, result_norepulsion.clusts))
# pointestimate_norepulsion = rcopy(R"""
#     library(salso)
#     salso($result_z_norepulsion)
#     """)

# Note - salso crashes in the above computation, so we use the MAP method for a point estimate instead. This should not really matter, since the coclustering matrix is basically singletons. 
pointestimate_norepulsion, _ = getpointestimate(result_norepulsion; method="MAP")
# Visualise the adjacency matrix of the point-estimate
sqmatrixplot(adjacencymatrix(pointestimate_norepulsion), 
title = "Adjacency Matrix of Point Estimate (No Repulsion)", colorbar=false)
savefig(joinpath(saveloc,"predicted_coclustering_matrix_norepulsion.pdf"))

# We can check the accuracy of the point estimate in terms of clustering metrics.
summarise(pointestimate_norepulsion, clusts)
# Save the results into a file
begin
local summary_file = open(joinpath(saveloc, "summary_norepulsion.txt"), "w")
show(summary_file, "text/plain", result_norepulsion)
println(summary_file, "")
summarise(summary_file, pointestimate_norepulsion, clusts)
close(summary_file)
end

# Next we try the Mixture of Finite Mixtures by J. Miller
# J. W. Miller and M. T. Harrison. Mixture models with a prior on the number of components. Journal of the American Statistical Association, Vol. 113, 2018, pp. 340-356.
# Julia package (BayesianMixtures) - https://github.com/jwmi/BayesianMixtures.jl
begin
    local tempmean = mean(points)
    points_MFM = [x .- tempmean for x in points]
    local tempsd = sqrt.(mean([x.*x for x in points_MFM]))
    points_MFM = [x ./ tempsd for x in points_MFM]
end
options_MFM = BayesianMixtures.options("MVNaaC", "MFM", points_MFM, numiters)
result_MFM = BayesianMixtures.run_sampler(options_MFM)

# Plot the posterior coclustering matrix
posterior_coclustering_MFM = zeros(Float64, length(points), length(points))
for i in 1:size(result_MFM.z, 2)
    posterior_coclustering_MFM .+= adjacencymatrix(Int.(result_MFM.z[:, i]))
end
posterior_coclustering_MFM ./= size(result_MFM.z, 2)
sqmatrixplot(posterior_coclustering_MFM,
title="Posterior Coclustering Probabilities (MFM)")
savefig(joinpath(saveloc,"MFM_posterior_coclustering_matrix.pdf"))

# Plot the posterior distribution of K:
histogram_pmf(result_MFM.t,
xlabel = "\$K\$", ylabel = "Probability", title = "Posterior Distribution of \$K\$ (MFM)")
savefig(joinpath(saveloc,"MFM_K_posterior.pdf"))

# Use the salso algorithm to find a point estimate
pointestimate_MFM = rcopy(R"""
    library(salso)
    salso($(permutedims(result_MFM.z)))
    """)
# Visualise the adjacency matrix of the point-estimate
sqmatrixplot(adjacencymatrix(pointestimate_MFM), 
    title = "Adjacency Matrix of Point Estimate (MFM)", colorbar=false)
savefig(joinpath(saveloc,"MFM_predicted_coclustering_matrix.pdf"))

# We can check the accuracy of the point estimate in terms of clustering metrics.
summarise(pointestimate_MFM, clusts)
# Save the results into a file
begin
    local summary_file = open(joinpath(saveloc, "summary_MFM.txt"), "w")
    summarise(summary_file, pointestimate_MFM, clusts)
    close(summary_file)
end

# Next we try a DPM model with normal kernel, also from the BayesianMixtures package
options_DPM = BayesianMixtures.options("MVNaaC", "DPM", points_MFM, numiters)
result_DPM = BayesianMixtures.run_sampler(options_DPM)

# Plot the posterior coclustering matrix
posterior_coclustering_DPM = zeros(Float64, length(points), length(points))
for i in 1:size(result_DPM.z, 2)
    posterior_coclustering_DPM .+= adjacencymatrix(Int.(result_DPM.z[:, i]))
end
posterior_coclustering_DPM ./= size(result_DPM.z, 2)
sqmatrixplot(posterior_coclustering_DPM,
title="Posterior Coclustering Probabilities (DPM)")
savefig(joinpath(saveloc,"DPM_posterior_coclustering_matrix.pdf"))

# Plot the posterior distribution of K:
histogram_pmf(result_DPM.t,
xlabel = "\$K\$", ylabel = "Probability", title = "Posterior Distribution of \$K\$ (DPM)")
savefig(joinpath(saveloc,"DPM_K_posterior.pdf"))

# Use the salso algorithm to find a point estimate
pointestimate_DPM= rcopy(R"""
    library(salso)
    salso($(permutedims(result_DPM.z)))
    """)
# Visualise the adjacency matrix of the point-estimate
sqmatrixplot(adjacencymatrix(pointestimate_DPM), 
    title = "Adjacency Matrix of Point Estimate (DPM)", colorbar=false)
savefig(joinpath(saveloc,"DPM_predicted_coclustering_matrix.pdf"))

# We can check the accuracy of the point estimate in terms of clustering metrics.
summarise(pointestimate_DPM, clusts)
# Save the results into a file
begin
    local summary_file = open(joinpath(saveloc, "summary_DPM.txt"), "w")
    summarise(summary_file, pointestimate_DPM, clusts)
    close(summary_file)
end

# Next we try k-means with the original points, setting K to the value obtained by the elbow method
pointestimate_kmeans = kmeans(makematrix(points), params.K_initial; maxiter=1000).assignments
sqmatrixplot(adjacencymatrix(pointestimate_kmeans), 
title = "Adjacency Matrix of Point Estimate (\$k\$-means)", colorbar=false)
savefig(joinpath(saveloc,"kmeans_predicted_coclustering_matrix.pdf"))
begin
    local summary_file = open(joinpath(saveloc, "summary_kmeans.txt"), "w")
    summarise(summary_file, pointestimate_kmeans, clusts)
    close(summary_file)
end

# Next we try k-medoids with the original points, setting K to the value obtained by the elbow method
pointestimate_kmedoids = kmedoids(distmatrix, params.K_initial; maxiter=1000).assignments
sqmatrixplot(adjacencymatrix(pointestimate_kmedoids), 
title = "Adjacency Matrix of Point Estimate (\$k\$-medoids)", colorbar=false)
savefig(joinpath(saveloc,"kmedoids_predicted_coclustering_matrix.pdf"))
begin
    local summary_file = open(joinpath(saveloc, "summary_kmedoids.txt"), "w")
    summarise(summary_file, pointestimate_kmedoids, clusts)
    close(summary_file)
end