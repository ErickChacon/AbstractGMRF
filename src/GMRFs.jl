module GMRFs

using Meshes
import SparseArrays: sparse, sparsevec, spdiagm, spzeros
import SuiteSparse # this is only for fixes
import FFTW
import Distributions: InverseGamma, Gamma, Distributions
import Random: AbstractRNG, randn!
import LinearAlgebra: cholesky, ldiv!, I, LinearAlgebra
import Graphs: Graphs, SimpleGraph, edges

include("utils.jl")
include("abstract-gmrf.jl")
# include("abstract-igmrf.jl")

export GMRF, RGMRF, CGMRF, GGMRF
export adjacency, structure, scale

end # module GMRFs
