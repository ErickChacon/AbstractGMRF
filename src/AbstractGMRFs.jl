module AbstractGMRFs

using Meshes
import SuiteSparse # this is only for fixes
import FFTW
import Distributions: InverseGamma, Gamma, Distributions
import Random: AbstractRNG, randn!
import LinearAlgebra: cholesky, ldiv!, I, LinearAlgebra

include("structure-matrix.jl")
include("gmrf.jl")
# include("igmrf.jl")

# export GMRF, RGMRF, CGMRF

end # module GMRF
