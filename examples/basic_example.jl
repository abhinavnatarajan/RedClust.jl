using RedClust
using Random: seed!
using StatsBase: counts
using CairoMakie

# Generate data 
seed!(44)
K = 10 # Number of clusters 
N = 100 # Number of points
data_σ = 0.25 # Variance of the normal kernel
data_dim = 10 # Data dimension
data = generatemixture(N, K; α = 10, σ = data_σ, dim = data_dim)
pnts, distM, clusts, probs, oracle_coclustering = data
# See the true adjacency matrix and oracle co-clustering matrices
fig1 = Figure()
Axis(fig1[1, 1], aspect = 1)
heatmap!(adjacencymatrix(clusts), colormap = :Blues)
fig1
fig2 = Figure()
Axis(fig2[1, 1], aspect = 1)
heatmap!(oracle_coclustering, colormap = :Blues)
fig2

# Determine the best the prior hyperparameters
params = fitprior(pnts, "k-means", false).params

# MCMC options
options = MCMCOptionsList()
data = MCMCData(D = distM)

# Run the sampler
result = runsampler(data, options, params)

# Check the results
# Posterior coclustering matrix
fig3 = Figure()
Axis(fig3[1, 1], aspect = 1)
heatmap!(result.posterior_coclustering, colormap = :Blues)
fig3
# Point-estimate adjacency matrix
fig4 = Figure()
Axis(fig4[1, 1], aspect = 1)
heatmap!(adjacencymatrix(result.pntestimate), colormap = :Blues)
fig4

# Posterior distribution of K
K_barplot = Figure(resolution = (400, 300), fontsize = 20)
ax = Axis(K_barplot[1, 1])
barplot!(ax, minimum(unique(result.K)):maximum(unique(result.K)), counts(result.K))
K_barplot
# Posterior distribution of r
r_hist = Figure(resolution = (400, 300), fontsize = 20)
ax = Axis(r_hist[1, 1])
hist!(ax, result.r, bins = 25)
r_hist
# Posterior distribution of p
p_hist = Figure(resolution = (400, 300), fontsize = 20)
ax = Axis(p_hist[1, 1])
hist!(ax, result.p, bins = 25)
p_hist
# Summary of MCMC and point estimate
summarise(result, clusts);

# Log-likelihood
loglik_plot = Figure(fontsize = 20)
ax = Axis(loglik_plot[1, 1])
lines!(ax, result.loglik)
loglik_plot
# Log-posterior
logposterior_plot = Figure(fontsize = 20)
ax = Axis(logposterior_plot[1, 1])
lines!(ax, result.logposterior)
logposterior_plot