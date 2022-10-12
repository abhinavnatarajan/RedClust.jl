using RedClust, Test
using StatsBase: sample

@testset "Convenience functions" begin

    @testset "makematrix" begin
        # check makematrix
        m = 500
        n = 1000
        temp = [rand(m) for i in 1:n]
        tempmatrix = makematrix(temp)
        @test sum([sum(tempmatrix[:, i] .== temp[i]) for i in 1:n]) == m * n
    end

    @testset "adjacencymatrix" begin
        # check adjacency matrix
        K = 20
        n = 500
        temp = sample(1:K, n)
        adjmatrix = adjacencymatrix(temp)
        tempsum = 0
        for (i, x) in pairs(adjmatrix)
            tempsum += (adjmatrix[i] == (temp[i[1]] == temp[i[2]]))
        end
        @test tempsum == n * n
    end

    @testset "sortlabels" begin
        # check sortlabels
        K = 20
        n = 500
        temp = sample(1:K, n)
        @test n^2 == sum(adjacencymatrix(temp) .== adjacencymatrix(sortlabels(temp)))
    end

end

K = 10 # Number of clusters 
N = 100 # Number of points
data_σ = 0.25 # Variance of the normal kernel
data_dim = 10 # Data dimension
data = generatemixture(N, K; α = 10, σ = data_σ, dim = data_dim)
pnts, distM, clusts, probs, oracle_coclustering = data

@testset begin 
    # check fitting the prior
    @testset "Fitting the prior" begin
        @test (
            try 
                fitprior(pnts, "k-means", false; useR = true) # kmeans with points, R
                true
            catch e
                false
            end
        )
        @test (
            try 
                fitprior(pnts, "k-medoids", false; useR = true) # kmedoids with points, R
                true
            catch e
                false
            end
        )
        @test (
            try 
                fitprior(distM, "k-medoids", true; useR = true) # kmedoids with distances, R
                true
            catch e
                false
            end
        )
        @test (
            try 
                fitprior(pnts, "k-means", false; useR = false) # kmeans with points, Julia
                true
            catch e
                false
            end
        )
        @test (
            try 
                fitprior(pnts, "k-medoids", false; useR = false) # kmedoids with points, Julia
                true
            catch e
                false
            end
        )
        @test (
            try 
                fitprior(distM, "k-medoids", true; useR = false) # kmedoids with distances, Julia
                true
            catch e
                false
            end
        )

        @test_throws ArgumentError fitprior(distM, "hierarchical", true; useR = false) # check that only 2 methods allowed

        @test_logs (:warn,) fitprior(pnts, "k-means", true; useR = false) # check that dist is disregarded when pnts is supplied
        
        @test_throws ArgumentError fitprior(distM, "k-means", true; useR = false) # check that kmeans does not work when distances are given

        @test_logs (:warn, ) fitprior(pnts, "k-means", false; Kmin = 0) # check that Kmin needs to be at least one

        @test_logs (:warn, ) fitprior(pnts, "k-means", false; Kmax = length(pnts) + 1) # check that Kmax needs to be at most N
    end


    # check the sampler 
    data = MCMCData(D = distM)
    @testset "Test sampler" begin
    @test (
        try 
            result = runsampler(data) # test defaults
            true
        catch e
            false
        end
    )
    end

    # Test point estimation
    result = runsampler(data)
    K = 50
    n = 500
    temp = sample(1:K, n)
    @testset "Test point estimation" begin
    @test infodist(temp, temp; normalised = true) ≈ 0 atol = 1e-9
    @test infodist(temp, temp; normalised = false) ≈ 0 atol = 1e-9
    @test binderloss(temp, temp; normalised = true) ≈ 0 atol = 1e-9
    @test binderloss(temp, temp; normalised = false) ≈ 0 atol = 1e-9
    @test (
        try 
            pntestimate1 = getpointestimate(result; method = "MAP")
            true
        catch e
            false
        end
    )
    @test (
        try 
            pntestimate2 = getpointestimate(result; method = "MLE")
            true
        catch e
            false
        end
    )
    @test (
        try 
            pntestimate3 = getpointestimate(result; loss = "binder", method = "MPEL")
            true
        catch e
            false
        end
    )
    @test (
        try 
            pntestimate4 = getpointestimate(result; loss = "omARI", method = "MPEL")
            true
        catch e
            false
        end
    )
    @test (
        try 
            pntestimate5 = getpointestimate(result; loss = "VI", method = "MPEL")
            true
        catch e
            false
        end
    )
    @test (
        try 
            pntestimate6 = getpointestimate(result; loss = "ID", method = "MPEL")
            true
        catch e
            false
        end
    )
    end
    
end