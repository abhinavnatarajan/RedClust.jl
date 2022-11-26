using RedClust
using Random: seed!
using StatsBase: counts
using Plots, StatsPlots
theme(:ggplot2)
default(fontfamily = "Computer Modern", 
guidefontsize = 16, 
tickfontsize = 16, 
colorbar_tickfontsize = 16, 
legend_font_pointsize = 16)

# Define convenience functions for plotting
function sqmatrixplot(X::Matrix, size_px = (500, 400))
    M, N = size(X)
    heatmap(
        X, 
        aspect_ratio=:equal, 
        c=:Blues, 
        xlim=(0,M), ylim=(0,N), 
        size = size_px)
end

function barplot(X::Vector)
    bar(sort(unique(X)), counts(X), legend = false, opacity = 0.8, linewidth = 0)
end

########## Generate data ##########
seed!(44)
K = 10 # Number of clusters 
N = 100 # Number of points
data_σ = 0.25 # Variance of the normal kernel
data_dim = 10 # Data dimension
data = generatemixture(N, K; α = 10, σ = data_σ, dim = data_dim);
points, distmatrix, clusts, probs, oracle_coclustering = data;
# Plot the true adjacency matrix and oracle co-clustering matrices as heatmaps
sqmatrixplot(adjacencymatrix(clusts))
sqmatrixplot(oracle_coclustering)
# Visualise distance matrix
sqmatrixplot(distmatrix)
# Histogram of distances
empirical_intracluster = uppertriangle(distmatrix)[
    uppertriangle(adjacencymatrix(clusts)) .== 1]
empirical_intercluster = uppertriangle(distmatrix)[
    uppertriangle(adjacencymatrix(clusts)) .== 0]
histogram(empirical_intercluster, opacity=0.7, linewidth=0, label="ICD", xlabel = "Distance", ylabel="Frequency")
histogram!(empirical_intracluster, opacity=0.7, bins=40, linewidth = 0, label="WCD")

############### MCMC ##############
# Determine the best the prior hyperparameters
params = fitprior(points, "k-means", false)
# Empirical vs prior predictive density of distances
pred_intracluster = sampledist(params, 10000, "intracluster")
pred_intercluster = sampledist(params, 10000, "intercluster")
density(pred_intracluster, label="Simulated WCD", xlabel = "Distance", ylabel = "Density", size = (700, 500), linewidth = 2, linestyle = :dash, linecolor = palette(:tab10)[1])
density!(empirical_intracluster, label="Empirical WCD", linewidth = 2, linecolor = palette(:tab10)[1])
density!(pred_intercluster, label="Simulated ICD", linewidth = 2, linestyle = :dash, linecolor = palette(:tab10)[2])
density!(empirical_intercluster, label="Empirical ICD", linewidth = 2, linecolor = palette(:tab10)[2])
# MCMC options
options = MCMCOptionsList(numiters=5000)
data = MCMCData(points)
# Run the sampler
result = runsampler(data, options, params)
pointestimate, index = getpointestimate(result; loss = "binder", method="MAP")

######## Check the results ########
# Summary of point estimate
summarise(pointestimate, clusts)
# Posterior coclustering matrix
sqmatrixplot(result.posterior_coclustering)
# Point-estimate adjacency matrix
sqmatrixplot(adjacencymatrix(pointestimate))
# Posterior distribution of K
bar(sort(unique(result.K)), counts(result.K), 
legend = false, 
opacity = 0.7, 
linewidth = 0, xlabel = "K", ylabel = "Frequency", 
size = (400, 400))
# Posterior distribution of r
histogram(result.r, opacity = 0.7, linewidth = 0, legend = false, ylabel = "Frequency", xlabel = "r")
# Posterior distribution of p
histogram(result.p, opacity = 0.7, linewidth = 0, legend = false, ylabel = "Frequency", xlabel = "p")
# Log-likelihood
plot(result.loglik, legend = false, xlabel = "Iteration", ylabel = "Log likelihood")
# Log-posterior
plot(result.logposterior, legend = false, xlabel = "Iteration", ylabel = "Log posterior")