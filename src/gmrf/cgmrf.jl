"""
    CGMRF(g::CartesianGrid, order::Integer, δ::Real, κ::Real)

Circulant Gaussian Markov random field over the regular grid `g` with scale parameter `κ` and structure
matrix defined by the penalty of n-th `order` and `δ` parameter.
"""
struct CGMRF <: AbstractGMRF
    g::CartesianGrid
    order::Integer
    δ::Number
    κ::Number
end

Base.length(d::CGMRF) = nelements(d.g)
scale(d::CGMRF) = d.κ
structure(d::CGMRF) = structure(d.g; δ = d.δ, order = d.order, circular = true)
structure_base(d::CGMRF) = structure_base(d.g; δ = d.δ, order = d.order)


## Random generator

function Distributions._rand!(rng::AbstractRNG, d::CGMRF, x::AbstractVector{T}) where T<:Real
    # get dimensions
    dims = size(d.g)
    n = prod(dims)

    # compute eigenvalues
    base = structure_base(d)
    λ = FFTW.fft(base)

    # simulate using fft
    z = randn(rng, reverse(dims)) + randn(rng, reverse(dims)) * im
    xaux = real(FFTW.fft(λ .^ (-1/2) .* z))
    if length(dims) == 2
        xaux = transpose(reverse(xaux, dims = 1))
    end
    copyto!(x, xaux)
    ldiv!(sqrt(n * scale(d)), x)

    return x
end

function Distributions._rand!(rng::AbstractRNG, d::CGMRF, x::AbstractArray{<:Real})
    dims = size(d.g)
    n = prod(dims)
    λ = FFTW.fft(structure_base(d))
    @inbounds for xi in Distributions.eachvariate(x, Distributions.variate_form(typeof(d)))
        z = randn(rng, reverse(dims)) + randn(rng, reverse(dims)) * im
        xaux = real(FFTW.fft(λ .^ (-1/2) .* z))
        if length(dims) == 2
            xaux = transpose(reverse(xaux, dims = 1))
        end
        copyto!(xi, xaux)
        ldiv!(sqrt(n * scale(d)), xi)
    end
    return x
end

## Logarithm of the pdf

function Distributions._logpdf(d::CGMRF, x::AbstractVector{<:Real})
    # get dimensions
    dims = size(d.g)
    n = prod(dims)

    # compute eigenvalues
    base = structure_base(d)
    λ = FFTW.fft(base)

    # compute u: vec(u) = S * vec(x)
    xmat = reverse(transpose(reshape(x, dims)), dims = 1)
    u = FFTW.fft(λ .* FFTW.ifft(xmat))

    # compute likelihood
    logpdf = -0.5 * n * log(2.0 * pi)
    logpdf += 0.5 * n * log(scale(d))
    logpdf += 0.5 * sum(log.(real(λ)))
    logpdf -= 0.5 * scale(d) * sum(xmat .* u)
    println("GMRF:")
    return real(logpdf)
end

# TODO: Find a way to not repeat this function from Distributions
@inline function Distributions.logpdf(d::CGMRF, x::AbstractArray{<:Real,1})
    @boundscheck begin
        size(x) == size(d) ||
            throw(DimensionMismatch("inconsistent array dimensions"))
    end
    return Distributions._logpdf(d, x)
end

# TODO: Custom modification from to Distributions.logpdf
@inline function Distributions.logpdf(d::CGMRF, x::AbstractArray{<:Real,M}) where {M}
    @boundscheck begin
        M > 1 ||
            throw(DimensionMismatch(
                "number of dimensions of `x` ($M) must be greater than number of dimensions of `d` (1)"
            ))
        ntuple(i -> size(x, i), Val(1)) == size(d) ||
            throw(DimensionMismatch("inconsistent array dimensions"))
    end

    dims = size(d.g)
    n = prod(dims)
    base = structure_base(d)
    λ = FFTW.fft(base)

    function auxfun(xi)
        xmat = reverse(transpose(reshape(xi, dims)), dims = 1)
        u = FFTW.fft(λ .* FFTW.ifft(xmat))
        sum(xmat .* u)
    end

    logpdf = -0.5 * n * log(2.0 * pi)
    logpdf += 0.5 * n * log(scale(d))
    logpdf += 0.5 * sum(log.(real(λ)))
    lpdf = @inbounds map(xi -> -0.5 * scale(d) * auxfun(x),
               Distributions.eachvariate(x, Distributions.variate_form(typeof(d))))
    println("GMRF:")
    return real(logpdf .+ lpdf)
end


