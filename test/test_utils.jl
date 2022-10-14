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

@testset "matsum and vecsum" begin
    x = rand(500, 500)
    v = rand(500)
    inds1 = sample(1:500, 200)
    inds2 = sample(1:500, 150)
    @test sum(x) ≈ RedClust.matsum(x)
    @test sum(x[inds1, inds2]) ≈ RedClust.matsum(x, inds1, inds2)
    @test sum(v) ≈ RedClust.vecsum(v)
    @test sum(v[inds1]) ≈ RedClust.vecsum(v, inds1)
end