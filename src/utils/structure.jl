"""
Return the structure matrix (S) of specified `order` associated to a `CartesianGrid`. The
structure matrix is defined such as Q = κS, where `Q` is the precision matrix of a GMRF or
IGMRF. An additional paramater δ can be defined such as S = D'D + δ*I, where D is a
difference matrix.
"""
function structure(g::CartesianGrid; δ = 0, order = 1, circular = false)
    D₁ = difference(g, circular = circular)
    if order == 1
        S = D₁'D₁
    elseif order > 1
        S = D₁'D₁ * structure(g, order = order - 1, circular = circular)
    end
    S + δ * I
end

function structure(g::SimpleGraph; δ = 0, order = 1)
    D₁ = difference(g, order = 1)
    if order == 1
        S = D₁'D₁
    elseif order > 1
        S = D₁'D₁ * structure(g, order = order - 1)
    end
    S + δ * I
end

"""
Return the base of a structure matrix (S) of specified `order` associated to a
`CartesianGrid`. The base of a structure matrix is defined such as Q = κS, where `Q` is the
circulant precision matrix of a GMRF or IGMRF. An additional paramater δ can be defined
such as S = D'D + δ*I, where D is a difference matrix.
"""
function structure_base(g::CartesianGrid{1}; δ = 0, order = 1)
    if order == 1
        base = sparsevec([1, 2, length(g)], [2.0, -1, -1])
    elseif order > 1
        w = centered([-1, 2, -1])
        base = imfilter(structure_base(g, order = order - 1), w, "circular")
    end
    base[1,1] += δ
    base
end

function structure_base(g::CartesianGrid{2}; δ = 0, order = 1)
    if order == 1
        (n1, n2) = size(g)
        base = sparse([1, 1, 2, 1, n2], [1, 2, 1, n1, 1], [4.0, -1, -1, -1, -1])
    elseif order > 1
        w = centered(sparse([0 -1 0; -1 4.0 -1; 0 -1 0]))
        base = imfilter(structure_base(g, order = order - 1), w, "circular")
    end
    base[1,1] += δ
    base
end
