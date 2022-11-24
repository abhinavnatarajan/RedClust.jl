using RedClust
using Random: seed!
using StatsBase: counts
using Plots
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
        fontfamily="Computer Modern",
        colorbar_tickfontsize=18,
        size = size_px)
end

function barplot(X::Vector)
    bar(sort(unique(X)), counts(X), legend = false)
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

############### MCMC ##############
# Determine the best the prior hyperparameters
params = fitprior(points, "k-means", false)
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
histogram(result.r, legend=false)
# Posterior distribution of p
histogram(result.p, legend=false)
# Log-likelihood
plot(result.loglik, legend = false)
# Log-posterior
plot(result.logposterior, legend = false)