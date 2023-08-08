module AbstractGMRFs

using Meshes
import SparseArrays: sparse, sparsevec, spdiagm, spzeros
import SuiteSparse # this is only for fixes
import FFTW
import Distributions: InverseGamma, Gamma, Distributions
import Random: AbstractRNG, randn!
import LinearAlgebra: cholesky, ldiv!, I, LinearAlgebra
import Graphs: Graphs, SimpleGraph, edges

include("adjacency.jl")
include("structure-matrix.jl")
include("gmrf.jl")
# include("igmrf.jl")

export GMRF, RGMRF, CGMRF, GGMRF
export adjacency, structure, scale

end # module GMRF
