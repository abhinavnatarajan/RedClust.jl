using Random: seed!
using HDF5: h5open, close, create_group

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
Return a read-only handle to a HDF5 file that contains example datasets from the main paper (see the introduction in the documentation for a persistent link to the paper). You must remember to close the file once you are done reading from it. This function is provided for reproducibility purposes only; it is recommended to read the datasets via the convenience function [`example_dataset`](@ref). 
"""
function example_datasets()
    h5open(datafilename, "r")
end

"""
    example_dataset(n::Int)

Returns a named tuple containing the dataset from the n``^\text{th}`` simulated example in the main paper. This dataset can also be generated using the following code in Julia v1.8.2:
    using RedClust
    using Random: seed!
    seed!(44)
    K = 10 # Number of clusters 
    N = 100 # Number of points
    data_σ = [0.25, 0.2, 0.18][n] # Variance of the normal kernel
    data_dim = [10, 50, 10][n] # Data dimension
    data = generatemixture(N, K; α = 10, σ = data_σ, dim = data_dim);

Note however that the above code may produce different results on different versions of Julia, even if they differ only by a minor version, because the random number generator in Julia is not guaranteed to be stable across versions. 

See also [`generate_mixture`](@ref). 
"""
function example_dataset(n::Int)
    if n ∉ [1, 2, 3]
        throw(ArgumentError("n must be 1, 2, or 3."))
    end
    datafile = example_datasets()
    egdata = datafile["example" * string(n)]
    x = read(egdata["points"])
    points = [x[:, i] for i in 1:last(size(x))]
    clusts = read(egdata["cluster_labels"])
    distmatrix = read(egdata["distance_matrix"])
    probs = read(egdata["cluster_weights"])
    oracle_coclustering = read(egdata["oracle_coclustering_probabilities"])
    close(datafile)
    data = (points = points, distmatrix = distmatrix, clusts = clusts, probs = probs, oracle_coclustering = oracle_coclustering)
    return data
end