"""
    CGMRF(base::AbstractMatrix, κ::Real)

Construct a Gaussian Markov random field (GMRF) with zero mean and block `circulant`
precision matrix with specified `base` and scale parameter `κ`.

    CGMRF(domain::CartesianGrid, order::Integer, κ::Real, δ::Real)

Construct a Gaussian Markov random field (GMRF) over specified `domain`.

The precision matrix the circulant GMRF is determined as `κS`, with `κ` serving as the
scaling parameter. The structure matrix `S` is computed using `S = D'D + δI`, where `D`
represents the difference matrix assuming circular neighboring topology depending of
selected `order`. The parameter `δ` is introduced to address positive definiteness issues
within the precision matrix. Notably, higher values of `δ` result in a faster correlation
decay, a higher `order` yields smoother processes, and a greater `κ` leads to reduced
process variance. This GMRF construction is applicable only to `CartesianGrid`s.
"""
struct CGMRF <: AbstractGMRF
    base::AbstractArray
    κ::Real
end

# Constructors

CGMRF(domain::CartesianGrid, order::Integer, κ::Real, δ::Real) =
    CGMRF(structure_base(domain; order = order, δ = δ), κ)

# Methods

Base.length(d::CGMRF) = length(d.base)
scale(d::CGMRF) = d.κ
structure_base(d::CGMRF) = d.base
# structure(d::CGMRF) = structure(d.g; δ = d.δ, order = d.order, circular = true)

## Random generator

function Distributions._rand!(rng::AbstractRNG, d::CGMRF, x::AbstractVector{T}) where T<:Real
    # get dimensions
    dims = reverse(size(structure_base(d)))
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
    dims = reverse(size(structure_base(d)))
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
    dims = reverse(size(structure_base(d)))
    n = prod(dims)

    # compute eigenvalues
    base = structure_base(d)
    λ = FFTW.fft(base)

    # compute u: vec(u) = S * vec(x)
    if length(dims) == 2
        xmat = reverse(transpose(reshape(x, dims)), dims = 1)
    else
        xmat = x
    end
    u = FFTW.fft(λ .* FFTW.ifft(xmat))

    # compute likelihood
    logpdf = -0.5 * n * log(2.0 * pi)
    logpdf += 0.5 * n * log(scale(d))
    logpdf += 0.5 * sum(log.(real(λ)))
    logpdf -= 0.5 * scale(d) * sum(xmat .* u)
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

    dims = reverse(size(structure_base(d)))
    n = prod(dims)
    base = structure_base(d)
    λ = FFTW.fft(base)

    function auxfun(xi)
        if length(dims) == 2
            xmat = reverse(transpose(reshape(xi, dims)), dims = 1)
        else
            xmat = xi
        end
        u = FFTW.fft(λ .* FFTW.ifft(xmat))
        sum(xmat .* u)
    end

    logpdf = -0.5 * n * log(2.0 * pi)
    logpdf += 0.5 * n * log(scale(d))
    logpdf += 0.5 * sum(log.(real(λ)))
    lpdf = @inbounds map(xi -> -0.5 * scale(d) * auxfun(xi),
               Distributions.eachvariate(x, Distributions.variate_form(typeof(d))))
    return real(logpdf .+ lpdf)
end


## IO

function Base.show(io::IO, ::MIME"text/plain", d::CGMRF)
    println(io, d)
    println(io, "  S base: ", summary(structure_base(d)))
    print(io, "  κ: ", scale(d))
end
