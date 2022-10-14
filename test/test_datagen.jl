K = 10 # Number of clusters 
N = 100 # Number of points
data_σ = 0.25 # Variance of the normal kernel
data_dim = 10 # Data dimension
@test_nothrow global data = generatemixture(N, K; α = 10, σ = data_σ, dim = data_dim)
pnts, distM, clusts, probs, oracle_coclustering = data
@test eltype(pnts) == Vector{Float64}
@test length(pnts[1]) == data_dim
@test length(unique(clusts)) == K
@test length(pnts) == N
@test size(distM, 1) == size(distM, 2)
@test size(distM, 1) == N
@test distM == distM'
@test sum(probs) ≈ 1
@test size(oracle_coclustering, 1) == size(oracle_coclustering, 2)
@test size(oracle_coclustering, 1) == N