using RedClust, Test
using StatsBase: sample

@testset "Convenience functions" begin

    @testset "makematrix" begin
        # check makematrix
        m = 500
        n = 1000
        for i = 1:50
            temp = [rand(m) for i in 1:n]
            tempmatrix = makematrix(temp)
            @test sum([sum(tempmatrix[:, i] .== temp[i]) for i in 1:n]) == m * n
        end
    end

    @testset "adjacencymatrix" begin
        # check adjacency matrix
        K = 20
        n = 500
        for i = 1:50
            temp = sample(1:K, n)
            adjmatrix = adjacencymatrix(temp)
            tempsum = 0
            for (i, x) in pairs(adjmatrix)
                tempsum += (adjmatrix[i] == (temp[i[1]] == temp[i[2]]))
            end
            @test tempsum == n * n
        end
    end

    @testset "sortlabels" begin
        # check sortlabels
        K = 20
        n = 500
        for i = 1:50
            temp = sample(1:K, n)
            @test n^2 == sum(adjacencymatrix(temp) .== adjacencymatrix(sortlabels(temp)))
        end
    end

end
# check fitting the prior

# check the sampler 