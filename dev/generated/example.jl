# Here we demonstrate the usage of RedClust through an example.

# We begin by setting up the necessary includes.

using RedClust, Plots, StatsPlots
using Random: seed!
using StatsBase: counts
using LinearAlgebra: triu, diagind


# We define some colors optimized for color-blind individuals based on [this article](https://www.nature.com/articles/nmeth.1618) by Bang Wong in Nature.
function wong_colors(alpha = 1.0)
    colors = [
        RGB(0/255, 114/255, 178/255), # blue
        RGB(230/255, 159/255, 0/255), # orange
        RGB(0/255, 158/255, 115/255), # green
        RGB(204/255, 121/255, 167/255), # reddish purple
        RGB(86/255, 180/255, 233/255), # sky blue
        RGB(213/255, 94/255, 0/255), # vermillion
        RGB(240/255, 228/255, 66/255), # yellow
    ]
    @. RGBA{Float32}(red(colors), green(colors), blue(colors), alpha)
end
# Set plotting defaults
default(fontfamily = "Computer Modern",
color_palette = wong_colors(0.8),
gridlinewidth = 1,
framestyle = :box,
linecolor = :match,
linewidth = 0.5,
guidefontsize = 14,
tickfontsize = 12,
colorbar_tickfontsize = 12,
legend_font_pointsize = 12,
plot_titlefontsize = 14
)
# seed the default RNG so that documentation remains stable
seed!(44)

# Next we define some convenience functions for plotting.

# Heatmap of square matrix
function sqmatrixplot(X::Matrix; kwargs...)
    M, N = size(X)
    heatmap(
        X,
        aspect_ratio=:equal,
        color=:Blues,
        xlim=(1,M), ylim=(1,N),
        yflip = true, xmirror=true;
        kwargs...)
end

# Histogram with integer bins
function histogram_pmf(X::AbstractVector{<:Integer}; binwidth::Int=1, kwargs...)
    xmin = minimum(X)
    xmax = maximum(X)
    binedges = xmin:binwidth:xmax
    c = counts(X)./length(X)
    bincounts = map(sum, [c[(i-xmin+1):minimum([i-xmin+1+binwidth-1, xmax-xmin+1])] for i in binedges])
    bar(binedges, bincounts,
    linewidth = 0,
    legend = false,
    xticks = binedges; kwargs...)
end

# Combine two symmetric square matrices together into the upper and lower triangle of a square matrix
function combine_sqmatrices(lower::Matrix, upper::Matrix, diagonal::String = "lower")
    if size(lower)[1] != size(lower)[2]
        throw(ArgumentError("Argument `lower` must be square, has dimensions $(size(lower))."))
    end
    if size(upper)[1] != size(upper)[2]
        throw(ArgumentError("Argument `upper` must be a square matrix, has dimensions $(size(upper))."))
    end
    if !all(size(lower) .== size(upper))
        throw(ArgumentError("Arguments `lower` and `upper` must have the same size."))
    end
    if !(eltype(lower) <: eltype(upper)) && !(eltype(upper) <: eltype(lower))
        throw(ArgumentError("Arguments must have compatible entries, got $(eltype(lower)) and $(eltype(upper))."))
    end
    if diagonal ∉ ["lower", "upper"]
        throw(ArgumentError("Keyword argument `diagonal` must be either \"lower\" or \"upper\"."))
    end
    result = copy(lower)
    temp = trues(size(lower))
    upper_idx = triu(temp, 1)
    diagonal_idx = diagind(temp)
    result[upper_idx] .= upper[upper_idx]
    result[diagonal_idx] .= ((diagonal == "lower") ? lower : upper)[diagonal_idx]
    return result
end

# ## Generating Data
# We can generate some example data using the function `generatemixture`.

begin
    K = 10 # Number of clusters
    N = 100 # Number of points
    data_σ = 0.25 # Variance of the normal kernel
    data_dim = 10 # Data dimension
    α = 10 # parameter for Dirichlet prior on cluster weights
    data = generatemixture(N, K;
    α = α, σ = data_σ, dim = data_dim)
    points, distmatrix, clusts, probs, oracle_coclustering = data
end

# Alternatively, the function `example_dataset` can be used to retrieve the datasets used in the original RedClust paper.

begin
    data = example_dataset(1)
    points, distmatrix, clusts, probs, oracle_coclustering = data
end

# We can visualise the true adjacency matrix of the observations with respect to the true clusters that they were drawn from, as well as the oracle coclustering matrix. The latter is the matrix of co-clustering probabilities of the observations conditioned upon the data generation process. This takes into account full information about the cluster weights (and how they are generated), the mixture kernels for each cluster, and the location and scale parameters for these kernels.

sqmatrixplot(combine_sqmatrices(oracle_coclustering, 1.0 * adjacencymatrix(clusts)), title = "Adjacency vs Oracle Co-clustering Probabilities \n(upper right and lower left triangle)\n")

# We can visualise the matrix of pairwise distances between the observations.

sqmatrixplot(distmatrix, title = "Matrix of Pairwise Distances")

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

# RedClust includes the function `fitprior` to heuristically choose prior hyperparameters based on the data.

params = fitprior(points, "k-means", false)

# We can check how good the chosen prior hyperparameters are by comparing the empirical distribution of distances to the (marginal) prior predictive distribution.

begin
    pred_intracluster = sampledist(params, "intracluster", 10000)
    pred_intercluster = sampledist(params, "intercluster", 10000)
    density(pred_intercluster,
    label="Simulated ICD", xlabel = "Distance", ylabel = "Density",
    linewidth = 2, linestyle = :dash)
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

# We can also evaluate the prior hyperparameters by checking the marginal predictive distribution on ``K`` (the number of clusters).

begin
    Ksamples = sampleK(params, 10000, N)
    histogram_pmf(Ksamples, legend = false,
    xticks=collect(0:10:maximum(Ksamples)),
    xlabel = "\$K\$", ylabel = "Probability", title = "Marginal Prior Predictive Distribution of \$K\$")
end

# Running the MCMC is straightforward. We set up the MCMC options using `MCMCOptionsList`.

options = MCMCOptionsList(numiters=50000)

# We then set up the input data using `MCMCData`.

data = MCMCData(points)

# We can then run the sampler using `runsampler`.

result = runsampler(data, options, params)

# The MCMC result contains several details about the MCMC, including acceptance rate, runtime, and convergence diagnostics. For full details see `MCMCResult`. In this example we have the ground truth cluster labels, so we can evaluate the result. For example, we can compare the posterior coclustering matrix to the oracle co-clustering probabilities.

sqmatrixplot(combine_sqmatrices(result.posterior_coclustering, oracle_coclustering),
title="Posterior vs Oracle Coclustering Probabilities")

# Plot the posterior distribution of ``K``:

histogram_pmf(result.K,
xlabel = "\$K\$", ylabel = "Probability", title = "Posterior Distribution of \$K\$")

# Plot the posterior distribution of ``r``:

begin
    histogram(result.r, normalize = :pdf,
    legend_font_pointsize=12,
    label="Empirical density", ylabel = "Density", xlabel = "\$r\$",
    title = "Posterior Distribution of \$r\$")
    density!(result.r,
    color=:black, linewidth = 2, linestyle=:dash,
    label="Kernel estimate", legend_font_pointsize=12)
end

# Plot the posterior distribution of ``p``:

begin
    histogram(result.p, normalize = :pdf,
    ylabel = "Density", xlabel = "\$p\$",
    title = "Posterior Distribution of \$p\$",
    label = "Empirical density",
    legend_font_pointsize=12)
    density!(result.p, color=:black, linewidth = 2, linestyle=:dash,
    label = "Kernel estimate",
    legend_position = :topleft)
end

# Plot the traceplot of the autocorrelation function of ``K``:

plot(result.K_acf, legend = false, linewidth = 1,
xlabel = "Lag", ylabel = "Autocorrelation",
title = "Autocorrelation Function of \$K\$")

# Plot the traceplot of the autocorrelation function of ``r``:

plot(result.r_acf, legend = false, linewidth = 1,
xlabel = "Lag", ylabel = "Autocorrelation",
title = "Autocorrelation Function of \$r\$")

# Plot the traceplot of the autocorrelation function of ``p``:

plot(result.p_acf, legend = false, linewidth = 1,
xlabel = "Lag", ylabel = "Autocorrelation",
title = "Autocorrelation Function of \$p\$")

# Check the trace plot of the log-likelihood to make sure the MCMC is moving well:

plot(result.loglik, legend = false, linewidth = 1,
xlabel = "Iteration", ylabel = "Log likelihood",
title = "Log-Likelihood Trace Plot")

# Check the trace plot of the log-posterior:

plot(result.logposterior, legend = false, linewidth = 1,
xlabel = "Iteration", ylabel = "Log posterior",
title = "Log-Posterior Trace Plot")

# The function `getpointestimate` finds an optimal point estimate, based on some notion of optimality. For example, to get the maximum a posteriori estimate we can run the following.

pointestimate, index = getpointestimate(result; method="MAP")

# We can compare the point-estimate to the true clustering through their adjacency matrices.

sqmatrixplot(combine_sqmatrices(adjacencymatrix(pointestimate), adjacencymatrix(clusts)),
title = "True Clustering vs MAP Point Estimate")

# We can check the accuracy of the point estimate in terms of clustering metrics.

summarise(pointestimate, clusts)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

