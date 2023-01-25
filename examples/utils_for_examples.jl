using LinearAlgebra: triu, diagind
# Set plotting defaults

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
function histogram_pmf(X::AbstractVector{<:Integer}; kwargs...)
    xvals = minimum(X):maximum(X)
    yvals = counts(X)./length(X)
    bar(xvals, yvals,
    linewidth = 0,
    legend = false,
    xticks = xvals; kwargs...)
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
