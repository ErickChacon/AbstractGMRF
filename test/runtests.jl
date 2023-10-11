using GMRFs
using Test
using Meshes
using SparseArrays
using LinearAlgebra
import Graphs
using Distributions

# Custom functions

xoxt(x::SparseMatrixCSC) = x .| x'
xpxt(x::SparseMatrixCSC) = x + x'

# Run tests

files = ["utils.jl", "gmrf.jl"]
# files = ["utils.jl"]

for file in files
    include(file)
end
