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
function sqmatrixplot(X::Matrix, size_px = (500, 400); kwargs...)
    M, N = size(X)
    heatmap(
        X, 
        aspect_ratio=:equal, 
        c=:Blues, 
        xlim=(0,M), ylim=(0,N), 
        size = size_px;
        kwargs...)
end

function integer_histogram(X; kwargs...)
    histogram(X, linewidth = 0, opacity = 0.7, legend = false, bins = minimum(X):maximum(X); kwargs...)
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
sqmatrixplot(adjacencymatrix(clusts), title="Adjacency Matrix")
sqmatrixplot(oracle_coclustering, title="Oracle Coclustering Probabilities")
# Visualise distance matrix
sqmatrixplot(distmatrix)
# Histogram of distances
empirical_intracluster = uppertriangle(distmatrix)[
    uppertriangle(adjacencymatrix(clusts)) .== 1]
empirical_intercluster = uppertriangle(distmatrix)[
    uppertriangle(adjacencymatrix(clusts)) .== 0]
histogram(empirical_intercluster, opacity=0.7, linewidth=0, label="ICD", xlabel = "Distance", ylabel="Frequency", title = "Observed distribution of distances")
histogram!(empirical_intracluster, opacity=0.7, bins=40, linewidth = 0, label="WCD")

############### Prior fitting ##############
# Determine the best the prior hyperparameters
params = fitprior(points, "k-means", false)
# Empirical vs prior predictive density of distances
pred_intracluster = sampledist(params, 10000, "intracluster")
pred_intercluster = sampledist(params, 10000, "intercluster")
density(pred_intracluster, label="Simulated WCD", xlabel = "Distance", ylabel = "Density", size = (700, 500), linewidth = 2, linestyle = :dash, linecolor = palette(:tab10)[1])
density!(empirical_intracluster, label="Empirical WCD", linewidth = 2, linecolor = palette(:tab10)[1])
density!(pred_intercluster, label="Simulated ICD", linewidth = 2, linestyle = :dash, linecolor = palette(:tab10)[2])
density!(empirical_intercluster, label="Empirical ICD", linewidth = 2, linecolor = palette(:tab10)[2])
# Visualise the marginal distribution on K
Ksamples = sampleK(params, 10000, N)
density(Ksamples, xlabel = "K", ylabel = "Density", linewidth = 2, legend = false, title = "Marginal Prior Predictive Density of K")

############### MCMC ##############
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
integer_histogram(result.K, xlabel = "K", ylabel = "Frequency", 
size = (600, 400), title = "Posterior Distribution of K")
# Posterior distribution of r
histogram(result.r, opacity = 0.7, linewidth = 0, legend = false, ylabel = "Frequency", xlabel = "r", title = "Posterior Distribution of r")
# Posterior distribution of p
histogram(result.p, opacity = 0.7, linewidth = 0, legend = false, ylabel = "Frequency", xlabel = "p", title = "Posterior Distribution of p")
# Log-likelihood
plot(result.loglik, legend = false, xlabel = "Iteration", ylabel = "Log likelihood", title = "Log-Likelihood Trace Plot")
# Log-posterior
plot(result.logposterior, legend = false, xlabel = "Iteration", ylabel = "Log posterior", title = "Log-Posterior Trace Plot")