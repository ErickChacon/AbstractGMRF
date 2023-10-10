using GMRFs
using Test
using Meshes
using SparseArrays
using LinearAlgebra
using Graphs

# custom functions

xoxt(x::SparseMatrixCSC) = x .| x'
xpxt(x::SparseMatrixCSC) = x + x'

# run tests

files = ["utils.jl"]

for file in files
    include(file)
end
