using RedClust
using Random: seed!
using StatsBase: counts
using Plots, StatsPlots
theme(:ggplot2)
default(fontfamily = "Computer Modern")

# Define convenience functions for plotting
function sqmatrixplot(X::Matrix, size_px = (500, 400))
    M, N = size(X)
    heatmap(
        X, 
        aspect_ratio=:equal, 
        c=:Blues, 
        xlim=(0,M), ylim=(0,N), 
        colorbar_tickfontsize=18,
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
histogram(empirical_intercluster, opacity=0.7, linewidth=0, label="ICD")
histogram!(empirical_intracluster, opacity=0.7, bins=40, linewidth = 0, label="WCD")

############### MCMC ##############
# Determine the best the prior hyperparameters
params = fitprior(points, "k-means", false)
# Empirical vs prior predictive density of distances
pred_intracluster = sampledist(params, 10000, "intracluster")
pred_intercluster = sampledist(params, 10000, "intercluster")
density(pred_intracluster, label="Simulated WCD")
density!(pred_intercluster, label="Simulated ICD")
density!(empirical_intracluster, label="Empirical WCD")
density!(empirical_intercluster, label="Empirical ICD")
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
barplot(result.K)
# Posterior distribution of r
histogram(result.r, opacity = 0.7, linewidth = 0, legend = false)
# Posterior distribution of p
histogram(result.p, opacity = 0.7, linewidth = 0, legend = false)
# Log-likelihood
plot(result.loglik, legend = false)
# Log-posterior
plot(result.logposterior, legend = false)