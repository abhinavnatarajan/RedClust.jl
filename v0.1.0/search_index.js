var documenterSearchIndex = {"docs":
[{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/#Main-functions","page":"Reference","title":"Main functions","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [RedClust]\nPages   = [\"Prior.jl\", \"MCMC.jl\"]\nOrder = [:function]","category":"page"},{"location":"reference/#RedClust.fitprior","page":"Reference","title":"RedClust.fitprior","text":"fitprior(data, algo; diss = false, Kmin = 1, Kmax = Int(floor(size(data)[end] / 2), useR = true)\n\nDetermines the best prior hyperparameters from the data. A notional clustering is obtained using k-means or k-medoids, and the distances are split into within-cluster distances and inter-cluster distances based on the notional clustering. These distances are then used to fit the prior hyperparameters using MLE and empirical Bayes sampling.   \n\nArguments\n\ndata::Union{Vector{Vector{Float64}}, Matrix{Float64}}: can either be a vector of (possibly multi-dimensional) observations, or a matrix with each column an observation, or a square matrix of pairwise dissimilarities. \nalgo::String: must be one of \"k-means\" or \"k-medoids\".\ndiss::bool = false: if true, data will be assumed to be a pairwise dissimilarity matrix. \nKmin::Integer: minimum number of clusters.\nKmax::Integer = Int(floor(size(data)[end] / 2)): maximum number of clusters. If left unspecified, it is set to half the number of observations.\nuseR::Bool = false: if false, will use the kmeans or kmedoids from the Julia package Clustering.jl. If true, will use kmeans or cluster::pam in R. \n\nReturns\n\nA named tuple containing the following\n\nparams::PriorHyperparamsList: the fitted hyperparameters.\nK::Integer: the number of clusters in the notional clustering.\nnotionalclustering::ClustLabelVector: the cluster labels in the notional clustering.\n\nSee also\n\nPriorHyperparamsList\n\n\n\n\n\n","category":"function"},{"location":"reference/#RedClust.runsampler","page":"Reference","title":"RedClust.runsampler","text":"runsampler(data[, options, params])\n\nRuns the MCMC sampler on the data.\n\nArguments\n\ndata::MCMCData: contains the distance matrix. \noptions = MCMCOptionsList(): contains the number of iterations, burnin, etc. \nparams = PriorHyperparamsList(): contains the prior hyperparameters for the model.\n\nReturns\n\nA struct of type MCMCResult containing the MCMC samples, convergence diagnostics, and summary statistics.\n\nSee also\n\nMCMCData, MCMCOptionsList, fitprior, PriorHyperparamsList.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Summary-and-clustering-evaluation","page":"Reference","title":"Summary and clustering evaluation","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [RedClust]\nPages   = [\"SummaryFunctions.jl\"]\nOrder = [:function]","category":"page"},{"location":"reference/#RedClust.evaluateclustering-Tuple{Vector{Int64}, Vector{Int64}}","page":"Reference","title":"RedClust.evaluateclustering","text":"evaluateclustering(clusts::ClustLabelVector, truth::ClustLabelVector)\n\nReturns a tuple of values quantifying the accuracy of the clustering assignments in clusts with respect to the ground truth clustering assignments in truth. \n\nReturn values\n\nbloss: Normalised Binder loss (= 1 - rand index). Lower is better. \nari: Adjusted Rand index. Higher is better. \nvidist: Variation of Information distances. Lower is better. \nvmeas: V-Measure. Higher is better. \nnmi: Normalised Mutual Information. Higher is better. \n\n\n\n\n\n","category":"method"},{"location":"reference/#RedClust.summarise","page":"Reference","title":"RedClust.summarise","text":"summarise(result::MCMCResult, truth::ClustLabelVector = [], printoutput = true)\n\nReturns a string that summarises the MCMC output. If truth is a vector of cluster labels, the output also contains clustering metrics to quantify the accuracy of the point-estimate from the MCMC with respect to truth. If printoutput = true, the output will be printed to console before being returned. \n\n\n\n\n\n","category":"function"},{"location":"reference/#Convenience-functions","page":"Reference","title":"Convenience functions","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [RedClust]\nPages   = [\"Utility.jl\"]\nOrder = [:function]","category":"page"},{"location":"reference/#RedClust.adjacencymatrix-Tuple{Vector{Int64}}","page":"Reference","title":"RedClust.adjacencymatrix","text":"adjacencymatrix(clusts::ClustLabelVector)\n\nReturns the nxn adjacency matrix corresponding to the given cluster label vector clusts, where n = length(clusts). \n\n\n\n\n\n","category":"method"},{"location":"reference/#RedClust.generatemixture-Tuple{Any, Any}","page":"Reference","title":"RedClust.generatemixture","text":"generatemixture(N, K; α = K, dim = K, radius = 1, σ = 0.1)\n\nGenerates a multivariate Normal mixture, with kernel weights generated from a Dirichlet prior. The kernels are centred at the vertices of a dim-dimensional simplex with edge length radius.\n\nArguments\n\nN::Integer: number of observations to generate.\nK::Integer: number of mixture kernels.\nα::Float64 = K: parameter for the Dirichlet prior.\ndim::Integer = K: dimension of the observations.\nradius::Float64 = 1: radius of the simplex whose vertices are the kernel means.\nσ::Float64 = 0.1: variance of each kernel.\n\nReturns\n\nNamed tuple containing the following fields-\n\npnts::Vector{Vector{Float64}}: a vector of N observations.\ndistM::Matrix{Float64}: an NxN matrix of pairwise Euclidean distances between the observations.\nclusts::ClustLabelVector: vector of N cluster assignments.\nprobs::Float64: vector of K cluster weights generated from the Dirichlet prior, used to generate the observations.\noracle_coclustering::Matrix{Float64}: NxN matrix of co-clustering probabilities, calculated assuming full knowledge of the cluster centres and cluster weights.\n\n\n\n\n\n","category":"method"},{"location":"reference/#RedClust.makematrix-Union{Tuple{Array{Vector{T}, 1}}, Tuple{T}} where T","page":"Reference","title":"RedClust.makematrix","text":"Convert a vector of vectors into a matrix, where each vector becomes a column in the matrix.\n\n\n\n\n\n","category":"method"},{"location":"reference/#RedClust.pointestimate-Tuple{Vector{Vector{Int64}}}","page":"Reference","title":"RedClust.pointestimate","text":"pointestimate(clusts::Vector{ClustLabelVector}; loss = \"VI\", usesalso = false)\n\nComputes a point estimate from a vector of samples of cluster allocations by searching for a minimiser of the posterior expectation of some loss function. \n\nArguments\n\nloss::String: must be either \"VI\" or \"binder\". Determines the loss function as either the Variation of Information distance or Binder loss. \n\nusesalso::Bool: if true, the SALSO algorithm is used. If false, the search space is restricted to the list of samples passed to the function.\n\n\n\n\n\n","category":"method"},{"location":"reference/#RedClust.sortlabels-Tuple{Vector{Int64}}","page":"Reference","title":"RedClust.sortlabels","text":"sortlabels(x::ClustLabelVector)\n\nReturns a cluster label vector y such that x and y have the same adjacency structure and labels in y occur in sorted ascending order.\n\n\n\n\n\n","category":"method"},{"location":"reference/#Types","page":"Reference","title":"Types","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [RedClust]\nOrder = [:type]","category":"page"},{"location":"reference/#RedClust.ClustLabelVector","page":"Reference","title":"RedClust.ClustLabelVector","text":"Alias for Vector{Int}\n\n\n\n\n\n","category":"type"},{"location":"reference/#RedClust.MCMCData","page":"Reference","title":"RedClust.MCMCData","text":"Contains the pairwise dissimilarities. Construct by calling MCMCData(D), where D is a square matrix of pairwise dissimilarities.\n\n\n\n\n\n","category":"type"},{"location":"reference/#RedClust.MCMCOptionsList","page":"Reference","title":"RedClust.MCMCOptionsList","text":"List of options for running the MCMC. \n\nConstructor:     MCMCOptionsList(; [numiters, burnin, thin, numGibbs, numMH, pointestimation, usesalso])\n\nConstructor arguments\n\nnumiters::Integer = 5000: number of iterations to run.\nburnin::Integer = floor(0.2 * numiters): number of iterations to discard as burn-in.\nthin::Integer = 1: will keep every thin samples.\nnumGibbs:Integer = 1: number of intermediate Gibbs scans in the split-merge step.\nnumMH:Integer = 1: number of split-merge steps per MCMC iteration.\npointestimation::Bool = true: whether to compute a point-estimate clustering from the posterior samples.\nusesalso::Bool = false: whether to use the SALSO algorithm to compute a point estimate.\n\n\n\n\n\n","category":"type"},{"location":"reference/#RedClust.MCMCResult","page":"Reference","title":"RedClust.MCMCResult","text":"Struct containing MCMC samples. \n\nFields\n\noptions::MCMCOptionsList: options passed to the sampler. \nparams::PriorHyperparamsList: prior hyperparameters used by the sampler. \nclusts::Vector{ClustLabelVector}: contains the clustering allocations. clusts[i] is a vector containing the clustering allocation from the ith sample. \npntestimate::ClustLabelVector: contains the point-estimate clustering allocation. If options.pointestimation == false then this is a vector of zeros. \nposterior_coclustering::Matrix{Float64}: the posterior coclustering matrix. \nK::Vector{Int}: posterior samples of K, i.e., the number of clusters. K[i] is the number of clusters in clusts[i].\nr::Vector{Float64}: posterior samples of the parameter r.\np::Vector{Float64}: posterior samples of the parameter p.\nK_mean, r_mean, p_mean: posterior mean of K, r, and p respectively. \nK_variance, r_variance, p_variance: posterior variance of K, r, and p respectively. \nK_acf::Vector{Float64}, r_acf::Vector{Float64}, p_acf::Vector{Float64}: autocorrelation function for K, r, and p respectively. \nK_iac, r_iac, and p_iac: integrated autocorrelation coefficient for K, r, and p respectively. \nK_ess::Float64, r_ess::Float64, and p_ess::Float64: effective sample size for K, r, and p respectively.\nloglik::Vector{Float64}: log-likelihood for each sample. \nsplitmerge_splits: Boolean vector indicating the iterations when a split proposal was used in the split-merge step. Has length numMH * numiters (see MCMCOptionsList).\nsplitmerge_acceptance_rate: acceptance rate of the split-merge proposals. \nr_acceptances: Boolean vector indicating the iterations (including burnin and the thinned out iterations) where the Metropolis-Hastings proposal for r was accepted. \nr_acceptance_rate:: Metropolis-Hastings acceptance rate for r.\nruntime: total runtime for all iterations.\nmean_iter_time: average time taken for each iteration. \n\n\n\n\n\n","category":"type"},{"location":"reference/#RedClust.PriorHyperparamsList","page":"Reference","title":"RedClust.PriorHyperparamsList","text":"Contains the prior hyperparameters for the model. Call the constructor with the line \n\nPriorHyperparamsList(; [kwargs])\n\nConstructor arguments\n\nδ1::Float64 = 1: the parameter delta_1.\nα::Float64 = 1: the shape parameter alpha in the prior for each lambda_k.\nβ::Float64 = 1: the rate parameter beta in the prior for each lambda_k.\nδ2::Float64 = 1: the parameter delta_2.\nζ::Float64 = 1: the shape parameter zeta in the prior for each theta_kt.\nγ::Float64 = 1: the rate parameter gamma in the prior for each theta_kt.\nη::Float64 = 1: the shape parameter eta in the prior for r.\nσ::Float64 = 1: the rate parameter sigma in the prior for r.\nu::Float64 = 1: the parameter u in the prior for p.\nv::Float64 = 1: the parameter v in the prior for p.\nproposalsd_r::Float64 = 1: standard deviation of the truncated Normal proposal when sampling r via a Metropolis-Hastings step.\nK_initial::Int = 1: initial number of clusters to start the MCMC.\nrepulsion::Bool = true: whether to use the repulsive component of the likelihood. \nmaxK::Int = 0: maxmimum number of clusters to allow. If set to 0 then we can get up to n clusters, where n is the number of observations.\n\nSee also\n\nfitprior\n\n\n\n\n\n","category":"type"},{"location":"","page":"Introduction","title":"Introduction","text":"CurrentModule = RedClust","category":"page"},{"location":"#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"RedClust is a Julia package for Bayesian clustering of high-dimensional Euclidean data using pairwise dissimilarity information instead of the raw observations. It uses an MCMC sampler to generate posterior samples from the space of all possible clustering structures on the data. ","category":"page"},{"location":"#Installation","page":"Introduction","title":"Installation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"The package can be installed by typing ]add RedClust into the Julia REPL or by the usual method:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"using Pkg\nPkg.add(\"RedClust\")","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"RedClust also requires R and the R package salso. If R is already installed, make sure the R_HOME environment variable is set to the R home directory (you could run R.home() in R to determine the location of this directory). If R or salso are not found, they are automatically installed during package installation.   ","category":"page"},{"location":"#Basic-example","page":"Introduction","title":"Basic example","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"using RedClust\n# Generate data\npnts, distM, clusts, probs, oracle_coclustering = \n\tgeneratemixture(N, K; α = 10, σ = data_σ, dim = data_dim)\n# Let RedClust choose the best prior hyperparameters\nparams = fitprior(pnts, \"k-means\", false).params\n# Set the MCMC options\noptions = MCMCOptionsList(numiters = 5000)\ndata = MCMCData(D = distM)\n# Run the sampler\nresult = runsampler(data, options, params)","category":"page"},{"location":"#Model","page":"Introduction","title":"Model","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"RedClust implements the model described in Natarajan et al. (2022). The key features are-","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The use of a random partition model with an unknown number of clusters K which allows for posterior inference on K. That is, there is a prior distribution on the space of all possible clustering structures with any number of clusters from one to the number of observations. The number of clusters is an object of inference to be determined by MCMC sampling. \nThe pairwise dissimilarities between observations are assumed to be mutually independent given the clustering assignment; that is, the clustering likelihood is a composite or pseduo-likelihood. \nThe clustering likelihood is comprised of a cohesive part and a repulsive part. ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The prior on the clustering structure can be chosen according to application in question. The current version of RedClust implements a microclustering prior (Miller et al., 2015; Betancourt et al., 2022). This means that the partition is generated by drawing cluster sizes n_1 ldots n_K from a random distribution nu (conditional upon n_1 + ldots n_K = n where n is the number of observations), and cluster labels are given by a uniform random permutation of ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"underbrace1 ldots 1_n_1 text times ldots underbraceK ldots K_n_K text times","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"RedClust implements a microclustering prior with nu a shifted negative binomial with random parameters r and p. ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"beginalign*\nr sim mathrmGamma(eta sigma)\np sim mathrmBeta(u v)\npi(rho_n mid r p) propto K p^n-K(1-p)^rKGamma(r)^-Kprod_k=1^K n_k Gamma(n_k+r-1)\nendalign*","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"where rho_n is the partition and K is the number of clusters in the partition. The partition likelihood is given by","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"beginalign*\nlambda_k oversetmathrmiidsim mathrmGamma(alpha beta) quad 1 leq k leq K\ntheta_kt oversetmathrmiidsim mathrmGamma(zeta gamma) quad 1 leq k  t leq K\nendalign*","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"pi(D mid rho_n boldsymbollambda boldsymboltheta r p) = left prod_k=1^K prod_substacki j in C_k  i neq j fracD_ij^delta_1 - 1lambda_k^delta_1Gamma(delta_1) exp(-lambda_k D_ij) right leftprod_substackt k = 1  t neq K prod_substacki in C_k  j in C_t fracD_ij^delta_2-1theta_kt^delta_2Gamma(delta_2) exp(-theta_ktD_ij) right","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"where boldsymbol D is the matrix of pairwise dissimilarities between observations.","category":"page"},{"location":"#Point-estimation","page":"Introduction","title":"Point estimation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"A clustering point-estimate boldsymbol c^* can be determined by searching for a clustering that minimises the expectation of the Variation of Information distance (Wade and Ghahramani, 2018) from a clustering boldsymbol c chosen randomly from the posterior distribution. That is, ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"boldsymbol c^* = argmin_boldsymbol c mathbbE_boldsymbol cd_VI(boldsymbol c boldsymbol c)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"A naive method is to restrict the search space to those clusterings visited by the MCMC sampler. \nA better method is the SALSO algorithm (Dahl et al., 2022), implemented in the R package salso. ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The MCMC sampler in RedClust automatically computes a point estimate using one of the methods above. The choice of method can be specified in the options passed to the sampler. See MCMCOptionsList.","category":"page"},{"location":"#Citing-this-package","page":"Introduction","title":"Citing this package","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"If you want to use this package in your work, please cite it as:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Natarajan, A., De Iorio, M., Heinecke, A., Mayer, E. and Glenn, S., 2022. ‘Cohesion and Repulsion in Bayesian Distance Clustering’, arXiv 2107.05414.","category":"page"},{"location":"#References","page":"Introduction","title":"References","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Betancourt, B., Zanella, G. and Steorts, R. C. (2022), ‘Random partition models for microclustering tasks’, Journal of the American Statistical Association 117(539), 1215–1227. DOI: 10.1080/01621459.2020.1841647.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Dahl, D. B., Johnson, D. J. and M¨uller, P. (2022), ‘Search algorithms and loss functions for bayesian clustering’, Journal of Computational and Graphical Statistics 0(0), 1–13. DOI: 10.1080/10618600.2022.2069779.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Miller, J., Betancourt, B., Zaidi, A., Wallach, H. and Steorts, R. C. (2015), ‘Microclustering: When the cluster sizes grow sublinearly with the size of the data set’, arXiv 1512.00792.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Wade, S. and Ghahramani, Z. (2018), ‘Bayesian Cluster Analysis: Point Estimation and Credible Balls (with Discussion)’, Bayesian Analysis 13(2), 559 – 626. DOI: 10.1214/17-BA1073.","category":"page"}]
}
