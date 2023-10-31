const datafilename = joinpath(@__DIR__, "..", "data", "example_datasets.h5")

# function generate_example_data(n::Int)
#     seed!(44)
#     K = 10 # Number of clusters
#     N = 100 # Number of points
#     data_σ = [0.25, 0.2, 0.18][n] # Variance of the normal kernel
#     data_dim = [10, 50, 10][n] # Data dimension
#     data = generatemixture(N, K; α = 10, σ = data_σ, dim = data_dim);
#     return data
# end

# function save_example_data()
#     datafile = h5open(datafilename, "w")
#     for i in 1:3
#         example = create_group(datafile, "example" * string(i))
#         data = generate_example_data(i)
#         points, distmatrix, clusts, probs, oracle_coclustering = data
#         example["points"] = makematrix(points)
#         example["distance_matrix"] = distmatrix
#         example["cluster_labels"] = clusts
#         example["cluster_weights"] = probs
#         example["oracle_coclustering_probabilities"] = oracle_coclustering
#     end
#     close(datafile)
# end

# save_example_data()

"""
Return a read-only handle to a HDF5 file that contains example datasets from the main paper. You must remember to close the file once you are done reading from it. This function is provided for reproducibility purposes only; it is recommended to read the datasets via the convenience function [`example_dataset`](@ref).
"""
function example_datasets()
    h5open(datafilename, "r")
end

@doc raw"""
	example_dataset(n::Int)

Returns a named tuple containing the dataset from the n``^{\mathrm{th}}`` simulated example in the main paper. This dataset was generated using the following code in Julia v1.8.1:
```julia
	using RedClust
	using Random: seed!
	seed!(44)
	K = 10 # Number of clusters
	N = 100 # Number of points
	data_σ = [0.25, 0.2, 0.18][n] # Variance of the normal kernel
	data_dim = [10, 50, 10][n] # Data dimension
	data = generatemixture(N, K; α = 10, σ = data_σ, dim = data_dim);
```

Note however that the above code may produce different results on your computer because the random number generator in Julia is not meant for reproducible results across different computers, different versions of Julia, or different versions of the Random.jl package, even with appropriate seeding. Therefore the datasets have been included with this package, and it is recommended to access them via this function.

See also [`generatemixture`](@ref).
"""
function example_dataset(n::Int)
    if n ∉ [1, 2, 3]
        throw(ArgumentError("n must be 1, 2, or 3."))
    end
    datafile = example_datasets()
    egdata = datafile["example"*string(n)]
    x = read(egdata["points"])
    points = [x[:, i] for i in 1:last(size(x))]
    clusts = read(egdata["cluster_labels"])
    distmatrix = read(egdata["distance_matrix"])
    probs = read(egdata["cluster_weights"])
    oracle_coclustering = read(egdata["oracle_coclustering_probabilities"])
    close(datafile)
    data = (points=points, distmatrix=distmatrix, clusts=clusts, probs=probs, oracle_coclustering=oracle_coclustering)
    return data
end
