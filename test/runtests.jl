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
    params = fitprior(pnts, "k-means", false; useR = false).params
    options1 = MCMCOptionsList(usesalso = true)
    options2 = MCMCOptionsList(usesalso = false)
    options3 = MCMCOptionsList(pointestimation = false, usesalso = false)
    data = MCMCData(D = distM)
    @testset "Test sampler" begin
    @test (
        try
            result1 = runsampler(data, options1, params) # point estimation with salso
            true
        catch e
            false
        end
    )
    @test (
        try
            result2 = runsampler(data, options2, params) # point estimation without salso 
            true
        catch e
            false
        end
    )
    @test (
        try
            result3 = runsampler(data, options3, params) # no point estimation
            sum(result3.pntestimate .== 0) == length(result3.pntestimate)
        catch e
            false
        end
    )
    @test (
        try 
            result4 = runsampler(data) # test defaults
            true
        catch e
            false
        end
    )
    end
end