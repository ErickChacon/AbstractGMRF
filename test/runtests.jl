using GMRFs
using Test
using Meshes
using SparseArrays
using LinearAlgebra
using Graphs

xxt(x::SparseMatrixCSC) = x .| x'

files = ["utils.jl"]

for file in files
    include(file)
end
