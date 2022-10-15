using RedClust, Test, Clustering
using StatsBase: sample

macro throws(expr)
    quote 
        try
            $(esc(expr))
            true
        catch
            false
        end
    end
end
macro test_nothrow(expr)
    quote 
        @test @throws $(esc(expr))
    end
end

@testset "master" begin
    @testset "utils" begin
        include("test_utils.jl")
    end
    @testset "datagen" begin
        include("test_datagen.jl")
    end
    @testset "fit prior" begin
        include("test_fitprior.jl")
    end
    @testset "sampler" begin
        include("test_sampler.jl")
    end
    @testset "point estimates" begin
        include("test_pointestimates.jl")
    end
end