module GMRFs

using Meshes
import SparseArrays: sparse, sparsevec, spdiagm, spzeros, findnz, SparseVector, SparseMatrixCSC
import SuiteSparse # this is only for fixes
import FFTW
import Distributions: InverseGamma, Gamma, Distributions
import Random: AbstractRNG, randn!
import LinearAlgebra: cholesky, ldiv!, I, LinearAlgebra
import Graphs: Graphs, SimpleGraph, edges, ne, adjacency_matrix
import ImageFiltering: imfilter, centered

include("utils.jl")
include("gmrf.jl")
# include("abstract-igmrf.jl")

export GMRF, RGMRF, CGMRF, GGMRF
export structure, scale

end # module GMRFs
