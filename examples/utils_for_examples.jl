using LinearAlgebra: triu, diagind
# Set plotting defaults
theme(:ggplot2)
default(fontfamily = "Computer Modern", 
guidefontsize = 16, 
tickfontsize = 16, 
colorbar_tickfontsize = 16, 
legend_font_pointsize = 16)

## Define convenience functions for plotting

# Heatmap of square matrix
function sqmatrixplot(X::Matrix; kwargs...)
    M, N = size(X)
    heatmap(
        X, 
        aspect_ratio=:equal, 
        c=:Blues, 
        xlim=(0,M), ylim=(0,N), 
        yflip = true, xmirror=true;
        kwargs...)
end

# Histogram with integer bins
function histogram_pmf(X::AbstractVector{<:Int}; kwargs...)
    bar(minimum(X):maximum(X), counts(X)./length(X), linewidth = 0, opacity = 0.7, legend = false; kwargs...)
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


