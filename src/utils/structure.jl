"""
Return the structure matrix (S) of specified `order` associated to a `CartesianGrid`. The
structure matrix is defined such as Q = κS, where `Q` is the precision matrix of a GMRF or
IGMRF. An additional paramater δ can be defined such as S = D'D + δ*I, where D is a
difference matrix.
"""
function structure(g::CartesianGrid; δ = 0, order = 1, cyclic = false)
    D = difference(g, order = order, cyclic = cyclic)
    D'D + δ * I
end

function structure(g::SimpleGraph; δ = 0, order = 1)
    D = difference(g, order = order)
    D'D + δ * I
end

"""
Return the base of a structure matrix (S) of specified `order` associated to a
`CartesianGrid. The base of a structure matrix is defined such as Q = κS, where `Q` is the
circulant precision matrix of a GMRF or IGMRF. An additional paramater δ can be defined
such as S = D'D + δ*I, where D is a difference matrix.
"""
function structure_base(g::CartesianGrid{1}; δ = 0, order = 1)
    n = g.dims[1]
    base = zeros(n)
    if order == 1
        base[[1, 2, end]] = [2 + δ, -1, -1]
    elseif order == 2
        base[[1, 2, 3, end-1, end]] = [6 + δ, -4, 1, 1, -4]
    else
        throw(ErrorException("not implemented"))
    end
    base
end

function structure_base(g::CartesianGrid{2}; δ = 0, order = 1)
    (n1, n2) = g.dims
    n = n1 * n2
    base = zeros(n2, n1)
    if order == 1
        base[1, 1] = 4 + δ
        base[1, 2] = -1
        base[2, 1] = -1
        base[1, end] = -1
        base[end, 1] = -1
    elseif order == 2
        base[1, 1] = 20 + δ
        base[1, 2] = -8
        base[2, 1] = -8
        base[1, end] = -8
        base[end, 1] = -8
        base[2, 2] = 2
        base[2, end] = 2
        base[end, 2] = 2
        base[end, end] = 2
        base[1, 3] = 1
        base[3, 1] = 1
        base[1, end-1] = 1
        base[end-1, 1] = 1
    else
        throw(ErrorException("not implemented"))
    end
    base
end
