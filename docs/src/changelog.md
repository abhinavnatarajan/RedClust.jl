# CHANGELOG

## [Unreleased]

### Added
- added log-posterior to result.
- added log-likelihood and log-posterior plots to basic example.
- point-estimation via maximum likelihood and maximum a posteriori.
- information distance between clusterings.
- add verbose option for [`runsampler`] and [`fitprior`].
- added bibliographic references to documentation.
- added changelog.

### Breaking
- pointestimate changed to `getpointestimate` with different signature.
- MCMCOptions constructor changed account for change in point-estimate calculation.  

### Fixed
- corrected computation of log-likelihood in result.

### Changed
- added bounds checking for `Kmin` and `Kmax` in [`fitprior`](@ref).
- added input message for better debugging in [`fitprior`](@ref).
- added progress bar for [`fitprior`](@ref) if `useR = false`. 
- added input validation for [`generatemixture`](@ref).
- added input validation for `binderloss` and `evaluateclustering`.
- separated computation of log-likelihood and log-prior.
- removed redundant loss function calculations when computing point-estimate.
- Binder loss calculation (`binderloss`) is now approximate (using randindex from Clustering.jl) but faster.
- [`runsampler`](@ref) no longer automatically calculates a point estimate. 
- [`summarise`](@ref) has a separate signature for MCMC output and for point estimate summary, and now provides more measures for point estimates.
- minor optimisations using `@inbounds` and `@simd`.

### To do

- remove redundant log-likelihood calculations in MCMC sampler.
- add options for using SALSO.
- add code coverage calculation for tests.