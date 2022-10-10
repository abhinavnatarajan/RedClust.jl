# CHANGELOG

## [Unreleased]

### Added
- added log-posterior to result.
- added log-likelihood and log-posterior plots to basic example.
- added changelog.

### Fixed
- added bounds checking for `Kmin` and `Kmax` in [`fitprior`](@ref).
- corrected computation of log-likelihood in result.


### Changed
- added input message for better debugging in [`fitprior`](@ref). 
- added progress bar for [`fitprior`](@ref) if `useR = false`. 
- separated computation of log-likelihood and log-prior.


### To do
- add verbose option for [`runsampler`] and [`fitprior`].
- remove redundant log-likelihood calculations in MCMC sampler.
- point-estimation via maximum likelihood and maximum a posteriori.
- add options for using SALSO.
- add code coverage calculation for tests.