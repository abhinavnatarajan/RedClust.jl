using RedClust, Plots, StatsPlots
using Random: seed!
using StatsBase: counts
include("utils_for_examples.jl")

## Load data
# Change example_dataset(1) to example_dataset(2) or example_dataset(3) for the other datasets
begin
    data = example_dataset(1)
    points, distmatrix, clusts, probs, oracle_coclustering = data;
    K = length(unique(clusts))
    N = length(points)
end

# Plot the true adjacency matrix and oracle co-clustering matrices as heatmaps
sqmatrixplot(adjacencymatrix(clusts), title = "Adjacency Matrix")
sqmatrixplot(oracle_coclustering, title = "Oracle Coclustering Probabilities")
# Visualise distance matrix
sqmatrixplot(distmatrix)
# Plot histogram of distances
begin 
    empirical_intracluster = uppertriangle(distmatrix)[
        uppertriangle(adjacencymatrix(clusts)) .== 1]
    empirical_intercluster = uppertriangle(distmatrix)[
        uppertriangle(adjacencymatrix(clusts)) .== 0]
    histogram(empirical_intercluster, 
    bins = minimum(empirical_intercluster):0.05:maximum(empirical_intercluster),
    label="ICD", xlabel = "Distance", ylabel="Frequency", 
    title = "Observed distribution of distances")
    histogram!(empirical_intracluster, 
    bins = minimum(empirical_intracluster):0.05:maximum(empirical_intracluster),
    label="WCD")
end

############### Prior fitting ##############
# Determine the best the prior hyperparameters
params = fitprior(points, "k-means", false)
# Plot Empirical vs prior predictive density of distances
begin 
    pred_intracluster = sampledist(params, "intracluster", 10000)
    pred_intercluster = sampledist(params, "intercluster", 10000)
    density(pred_intracluster, 
    label="Simulated WCD", xlabel = "Distance", ylabel = "Density", 
    size = (700, 500), 
    linewidth = 2, linestyle = :dash)
    density!(empirical_intracluster, 
    label="Empirical WCD", 
    linewidth = 2, primary = false)
    density!(pred_intercluster, 
    label="Simulated ICD", 
    linewidth = 2, linestyle = :dash)
    density!(empirical_intercluster, 
    label="Empirical ICD", 
    linewidth = 2, primary = false)
end
# Visualise the marginal distribution on K
begin
    Ksamples = sampleK(params, 10000, N)
    density(Ksamples, xlabel = "K", ylabel = "Density", linewidth = 2, legend = false, title = "Marginal Prior Predictive Density of K")
end


############### MCMC ##############
# MCMC options
begin
    options = MCMCOptionsList(
        numiters=5000, numGibbs = 50)
    data = MCMCData(points)
end

# Run the sampler
begin
    result = runsampler(data, options, params)
    pointestimate, index = getpointestimate(result; loss = "binder", method="MAP")
end

######## Check the results ########
# Summary of point estimate
summarise(pointestimate, clusts)
# Posterior coclustering matrix
sqmatrixplot(combine_sqmatrices(result.posterior_coclustering, oracle_coclustering), 
title="Posterior vs Oracle Coclustering Probabilities")
# Point-estimate adjacency matrix
sqmatrixplot(combine_sqmatrices(adjacencymatrix(pointestimate), adjacencymatrix(clusts)), 
title = "True Clustering vs MAP Point Estimate")
# Posterior distribution of K
histogram_pmf(result.K, xlabel = "K", ylabel = "PMF", 
size = (400, 400), title = "Posterior Distribution of K")
# Posterior distribution of r
begin
    histogram(result.r, normalize = :pdf,
    legend_font_pointsize=12, 
    label="Empirical density", ylabel = "Density", xlabel = "r", 
    title = "Posterior Distribution of r")
    density!(result.r, 
    color=:black, linewidth = 2, linestyle=:dash, 
    label="Kernel estimate", legend_font_pointsize=12)
end
# Posterior distribution of p
begin
    histogram(result.p, normalize = :pdf, 
    ylabel = "Density", xlabel = "p", 
    title = "Posterior Distribution of p", 
    label = "Empirical density", legend_font_pointsize=12)
    density!(result.p, color=:black, linewidth = 2, linestyle=:dash, 
    label = "Kernel estimate")
end
# Log-likelihood
plot(result.loglik, legend = false, linewidth = 1,
xlabel = "Iteration", ylabel = "Log likelihood", 
title = "Log-Likelihood Trace Plot")
# Log-posterior
plot(result.logposterior, legend = false, linewidth = 1,
xlabel = "Iteration", ylabel = "Log posterior", 
title = "Log-Posterior Trace Plot")