using CairoMakie

function sqmatrixplot(M::Matrix, size_inches = (4.5, 4), fontsize = 18)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = fontsize)
    Axis(fig[1, 1], aspect = 1)
    heatmap!(M, colormap = :Blues)
    return fig
end

function lineplot(v::Vector{Float64}, size_inches = (4, 3), fontsize = 18)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = fontsize)
    Axis(fig[1, 1])
    lines!(v)
    return fig
end