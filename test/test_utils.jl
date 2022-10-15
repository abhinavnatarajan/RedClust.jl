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

@testset "printing time" begin
    @test RedClust.prettytime(1e-9) == "1.00 ns"
    @test RedClust.prettytime(999e-9) == "999.00 ns"
    @test RedClust.prettytime(1e-6) == "1.00 μs"
    @test RedClust.prettytime(999e-6) == "999.00 μs"
    @test RedClust.prettytime(1e-3) == "1.00 ms"
    @test RedClust.prettytime(999e-3) == "999.00 ms"
    @test RedClust.prettytime(1) == "1.00 s"
    @test RedClust.prettytime(5) == "5.00 s"
    @test RedClust.prettytime(60) == "1 min"
    @test RedClust.prettytime(120) == "2 mins"
    @test RedClust.prettytime(65) == "1 min 5 s"
    @test RedClust.prettytime(125) == "2 mins 5 s"
    @test RedClust.prettytime(3600) == "1 hr"
    @test RedClust.prettytime(7200) == "2 hrs"
    @test RedClust.prettytime(7205) == "2 hrs 5 s"
    @test RedClust.prettytime(7265) == "2 hrs 1 min 5 s"
    @test RedClust.prettytime(7325) == "2 hrs 2 mins 5 s"
    @test RedClust.prettytime(24 * 3600) == "1 day"
    @test RedClust.prettytime(2 * 24 * 3600) == "2 days"
end