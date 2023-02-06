using Pkg
Pkg.activate(joinpath(@__DIR__))

# replace with the home directory of your R installation
ENV["R_HOME"] = "C:/Program Files/R/R-4.2.2/"

Pkg.add(url="https://github.com/jwmi/BayesianMixtures.jl") # this is an unregistered dependency, so we include it here until Julia is updated to allow such dependencies
Pkg.instantiate()
using RCall: rcopy, @R_str
salso_is_installed = rcopy(R"""
    require(salso)
    """)
salso_is_installed || error("The R package salso could not be loaded.")

@info "Successfully set up environment!"

