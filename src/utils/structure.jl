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
    n = size(g)
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
    if order == 1
        (n1, n2) = size(g)
        base = sparse([1, 1, 2, 1, n2], [1, 2, 1, n1, 1], [4.0, -1, -1, -1, -1])
    elseif order > 1
        w = centered(sparse([0 -1 0; -1 4.0 -1; 0 -1 0]))
        base = imfilter(structure_base(g; δ = 0, order = order - 1), w, "circular")
    end
    base[1,1] += δ
    base
end
