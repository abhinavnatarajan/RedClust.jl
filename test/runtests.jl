using RedClust, Test
using StatsBase: sample

@testset "Convenience functions" begin

# check makematrix
m = 500
n = 1000
temp = [rand(m) for i in 1:n]
tempmatrix = makematrix(temp)
@test sum([sum(tempmatrix[:, i] .== temp[i]) for i in 1:n]) == m * n

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

# check sortlabels
K = 20
n = 500
temp = sample(1:K, n)
@test n^2 == sum(adjacencymatrix(temp) .== adjacencymatrix(sortlabels(temp)))

end
# check fitting the prior

# check the sampler 