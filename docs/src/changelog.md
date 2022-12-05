# Changelog

## [1.0.0]

### Added
- examples from the main paper are now included in the examples folder. Data from the original paper is included in the data folder in the standard HDF5 format. 
- the examples now have more plots.
- functionality to sample the distances from the prior predictive distribution (see [`sampledist`](@ref)).
- functionality to sample ``K`` from its marginal prior predictive (see [`sampleK`](@ref)).
- [`generatemixture`](@ref) accepts a random number generator or a seed for the default RNG for reproducibility of results. 

### Removed 
- dependency on RCall and the various calls to the salso algorithm were removed. It is left to the user to make these calls if necessary. 

## [0.2.2]

### Fixed
- corrected time not printing in the correct format. 
- moved `prettytime` and `prettynumber` to `utils.jl`.

## [0.2.1]

### Added
- added tests for pretty printing.
- added tests for custom loss function in [`getpointestimate`](@ref). 

### Fixed
- verbose output in sampler missing newline. 
- custom loss function when computing a point estimate. 

### Removed
- `_infodist` not in use, removed.
  
## [0.2.0]

### Added
- added log-posterior to result.
- added log-likelihood and log-posterior plots to basic example.
- point-estimation via maximum likelihood and maximum a posteriori.
- [`infodist`](@ref) function to compute the information distance between clusterings.
- convenience constructor for [`MCMCData`](@ref).
- add verbose option for [`runsampler`](@ref) and [`fitprior`](@ref).
- pretty printing for [`MCMCData`](@ref), [`MCMCOptionsList`](@ref), and [`PriorHyperparamsList`](@ref).
- added bibliographic references to documentation.
- added changelog.

### Breaking
- MCMCOptions constructor changed account for change in point-estimate calculation. 
- [`fitprior`](@ref) now only returns the hyperparameter list.
- removed `summarise` for [`MCMCResult`](@ref) objects (use pretty printing instead).

### Fixed
- corrected computation of log-likelihood in result.
- fixed edge case error in [`fitprior`](@ref) when the found value of `K` is either 1 or `N`. This case arises only for edge values of `Kmin` and `Kmax`, since the elbow method will otherwise not choose 1 or `N` for `K`.

### Changed
- added bounds checking for `Kmin` and `Kmax` in [`fitprior`](@ref).
- added input message for better debugging in [`fitprior`](@ref).
- added progress bar for [`fitprior`](@ref) if `useR = false`. 
- added input validation for [`generatemixture`](@ref).
- added input validation for [`binderloss`](@ref) and [`evaluateclustering`](@ref).
- added input validation for the constructors of [`MCMCData`](@ref) and [`MCMCOptionsList`](@ref).
- separated computation of log-likelihood and log-prior.
- removed redundant loss function calculations when computing point-estimate.
- Binder loss calculation (`binderloss`) is now approximate (using randindex from Clustering.jl) but faster.
- [`runsampler`](@ref) no longer automatically calculates a point estimate. 
- [`summarise`](@ref) has a separate signature for MCMC output and for point estimate summary, and now provides more measures for point estimates.
- minor optimisations using `@inbounds`, `@simd`, and `@turbo`.
- using StaticArrays.jl and separate function for restricted Gibbs scan. 
- speedup via non-generic implementation of matrix and vector sums with LoopVectorization.jl.