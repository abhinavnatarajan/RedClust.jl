# CHANGELOG

## Dev

<!-- ### New features -->

### Bug fixes
- added bounds checking for `Kmin` and `Kmax` in [`fitprior`](@ref).

### Minor changes
- added input message for better debugging in [`fitprior`](@ref). 
- added progress bar for [`fitprior`](@ref) if `useR = false`. 

### Documentation 
- added changelog

### To do
- remove redundant log-likelihood calculations in MCMC sampler
- add log-posterior computation
- point-estimation via maximum likelihood and maximum a posteriori
- add options for using SALSO
- add code-coverage