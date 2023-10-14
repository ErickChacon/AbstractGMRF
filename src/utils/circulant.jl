"""
    spcirculant(x)

Return a sparse circulant or block-circulant matrix based on the input sparse vector or
matrix `x`.

The size of the resulting matrix is `n×n`, where `n` represents the total number of
elements in `x`. If `x` is a `SparseVector`, the function generates a `n×n` sparse
circulant matrix. Conversely, if `x` is a `SparseMatrixCSC`, the function creates a `n×n`
sparse block circulant matrix.
"""

function spcirculant(x::SparseVector)
    # get indices and values
    n = length(x)
    (is, vals) = findnz(x)

    # compute the indices for the full structure matrix
    Is = repeat([i for i in 1:n], length(vals))
    Js = [mod1(i -is[k]+1, n) for i in 1:n, k in 1:length(vals)]
    Vs = repeat(vals, inner = n)

    # spare structure matrix
    sparse(Is, vec(Js), Vs, n, n)
end

function spcirculant(x::SparseMatrixCSC)
    # get base and values
    (n1,n2) = size(x)
    (is, js, vals) = findnz(x)

    # compute the indices for the full structure matrix
    Is = repeat(vec([i + n1 * (b-1) for i in 1:n1, b in 1:n2]), length(vals))
    Js = [mod1(n1*(js[k]-1) + n1 * (b-1) + mod1(i -is[k]+1, n1), n1*n2) for i in 1:n1, b in 1:n2, k in 1:length(vals)]
    Vs = repeat(vals, inner = n1 * n2)

    # spare structure matrix
    sparse(Is, vec(Js), Vs, n1 * n2, n1 * n2)
end

