# Here we demonstrate the usage of RedClust through an example.
# We begin by setting up the necessary includes.
using Pkg
Pkg.activate(joinpath(@__DIR__))
Pkg.develop(Pkg.PackageSpec(path = joinpath(@__DIR__, "..")))
Pkg.instantiate()
using RedClust, Plots, StatsPlots
using Random: seed!
using StatsBase: counts
using LinearAlgebra: triu, diagind
include("utils_for_examples.jl")

dataset_number = 1

# seed the default RNG so that documentation remains stable
seed!(44)

# Folder to save things to
saveloc = joinpath(@__DIR__, "figures", "simulation_study_" * string(dataset_number))
if !isdir(saveloc)
    mkpath(saveloc)
end

# Get the data
begin
    data = example_dataset(dataset_number)
    points, distmatrix, clusts, probs, oracle_coclustering = data
    N = length(points)
end

# We can visualise the true adjacency matrix of the observations with respect to the true clusters that theWCDy were drawn from, as well as the oracle coclustering matrix. The latter is the matrix of co-clustering probabilities of the observations conditioned upon the data generation process. This takes into account full information about the cluster weights (and how they are generated), the mixture kernels for each cluster, and the location and scale parameters for these kernels.
sqmatrixplot(combine_sqmatrices(oracle_coclustering, 1.0 * adjacencymatrix(clusts)), title = "Adjacency vs Oracle Co-clustering Probabilities \n(upper right and lower left triangle)\n")
savefig(joinpath(saveloc,"adjacency_and_oracle_matrix.pdf"))

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
    bins = minimum(empirical_intercluster):0.05:maximum(empirical_intercluster),
    label="ICD", xlabel = "Distance", ylabel="Frequency",
    title = "Observed distribution of distances")
    histogram!(empirical_intracluster,
    bins = minimum(empirical_intracluster):0.05:maximum(empirical_intracluster),
    label="WCD")
end
savefig(joinpath(saveloc,"histogram_distances.pdf"))

# RedClust includes the function `fitprior` to heuristically choose prior hyperparameters based on the data.
params = fitprior(points, "k-means", false)

# We can check how good the chosen prior hyperparameters are by comparing the empirical distribution of distances to the (marginal) prior predictive distribution.
begin
    pred_intracluster = sampledist(params, "intracluster", 10000)
    pred_intercluster = sampledist(params, "intercluster", 10000)
    density(pred_intracluster,
    label="Simulated WCD", xlabel = "Distance", ylabel = "Density",
    linewidth = 2, linestyle = :dash)
    density!(empirical_intracluster,
    label="Empirical WCD",
    color = 1, linewidth = 2)
    density!(pred_intercluster,
    label="Simulated ICD",
    linewidth = 2, linestyle = :dash, color = 2)
    density!(empirical_intercluster,
    label="Empirical ICD",
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
options = MCMCOptionsList(numiters=50000)

# We then set up the input data using `MCMCData`.
data = MCMCData(points)

# We can then run the sampler using `runsampler`.
result = runsampler(data, options, params)

# The MCMC result contains several details about the MCMC, including acceptance rate, runtime, and convergence diagnostics. For full details see `MCMCResult`. In this example we have the ground truth cluster labels, so we can evaluate the result. For example, we can compare the posterior coclustering matrix to the oracle co-clustering probabilities.
sqmatrixplot(combine_sqmatrices(result.posterior_coclustering, oracle_coclustering),
title="Posterior vs Oracle Coclustering Probabilities")
savefig(joinpath(saveloc,"posterior_vs_oracle_coclustering.pdf"))

# Plot the posterior distribution of K:
histogram_pmf(result.K,
xlabel = "\$K\$", ylabel = "Probability", title = "Posterior Distribution of \$K\$")
savefig(joinpath(saveloc,"K_posterior.pdf"))

# Plot the posterior distribution of r:
begin
    histogram(result.r, normalize = :pdf,
    legend_font_pointsize=12, 
    label="Empirical density", ylabel = "Density", xlabel = "\$r\$", 
    title = "Posterior Distribution of \$r\$")
    density!(result.r, 
    color=:black, linewidth = 2, linestyle=:dash, 
    label="Kernel estimate", legend_font_pointsize=12)
end
savefig(joinpath(saveloc,"r_posterior.pdf"))

# Plot the posterior distribution of p:
begin
    histogram(result.p, normalize = :pdf,
    ylabel = "Density", xlabel = "\$p\$",
    title = "Posterior Distribution of \$p\$",
    label = "Empirical density",
    legend_font_pointsize=12,
    legend_position = :topleft)
    density!(result.p, color=:black, linewidth = 2, linestyle=:dash,
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

# The function `getpointestimate` finds an optimal point estimate, based on some notion of optimality. For example, to get the maximum a posteriori estimate we can run the following.
pointestimate, index = getpointestimate(result; method="MAP")

# We can compare the point-estimate to the true clustering through their adjacency matrices.
sqmatrixplot(combine_sqmatrices(adjacencymatrix(pointestimate), adjacencymatrix(clusts)),
title = "True Clustering vs MAP Point Estimate")
savefig(joinpath(saveloc,"predicted_vs_true_coclustering_matrix.pdf"))

# We can check the accuracy of the point estimate in terms of clustering metrics.
summarise(pointestimate, clusts)

# Next we try the Mixture of Finite Mixtures by J. Miller
# J. W. Miller and M. T. Harrison. Mixture models with a prior on the number of components. Journal of the American Statistical Association, Vol. 113, 2018, pp. 340-356.
# Julia package (BayesianMixtures) - https://github.com/jwmi/BayesianMixtures.jl

