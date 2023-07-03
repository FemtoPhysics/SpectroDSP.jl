"""
Cooley-Tukey butterfly computation.

    ctb!(ya::VecI{Complex{T}},
         xa::VecI{Complex{T}},
         wa::VecI{Complex{T}},
         si::Int, hs::Int, ns::Int,
         ss::Int, pd::Int) where T<:AbstractFloat
-----------------------------------
1. `ya`: destination array
2. `xa`: source array
3. `wa`: twiddle factors array
4. `si`: starting index
5. `hs`: half of FFT size
6. `ns`: number of steps
7. `ss`: step size of the current subproblem
8. `pd`: index difference of butterfly computation pair
"""
function ctb!(
        ya::VecI{Complex{T}},
        xa::VecI{Complex{T}},
        wa::VecI{Complex{T}},
        si::Int,
        hs::Int,
        ns::Int,
        ss::Int,
        pd::Int
    ) where T<:AbstractFloat
    wi = 1
    yi = xi = si
    for _ in 1:ns
        xj = xi + hs
        @inbounds begin
            ya[yi]    = (xa[xi] + xa[xj])
            ya[yi+pd] = (xa[xi] - xa[xj]) * wa[wi]
        end
        yi += ss
        xi += pd
        wi += pd
    end
    return nothing
end

"""
Decimation-in-time FFT with a naturally ordered input-output.

    ditnn!(sa::VecI{Complex{T}},
           ba::VecI{Complex{T}},
           wa::VecI{Complex{T}},
           hs::Int) where T<:AbstractFloat
-------------------------------------------------------------
1. `sa`: signal array
2. `ba`: buffer array
3. `wa`: twiddle factors array
4. `sf`: switch flag
"""
function ditnn!(sa::VecI{Complex{T}}, ba::VecI{Complex{T}},
                wa::VecI{Complex{T}}, hs::Int) where T<:AbstractFloat
    ns = hs
    pd = 1
    ss = 2
    sf = false

    while ns > 0
        if sf
            for si in 1:pd
                ctb!(sa, ba, wa, si, hs, ns, ss, pd)
            end
        else
            for si in 1:pd
                ctb!(ba, sa, wa, si, hs, ns, ss, pd)
            end
        end

        ns >>= 1
        pd <<= 1
        ss <<= 1
        sf = !sf
    end

    return nothing
end

"""
Radix2 Kernel for 1D FFT

Initialize and return a kernel with `fftsize` and specific type `T<:AbstractFloat`.

    Radix2FFT{T}(fftsize::Integer) where T<:AbstractFloat

Initialize and return a kernel with `fftsize` and default type `Float64`.

    Radix2FFT(fftsize::Integer)
"""
struct Radix2FFT{T<:AbstractFloat}
    cache  ::Vector{Complex{T}}
    twiddle::Vector{Complex{T}}
    fftsize::Int
    ifswap ::Bool

    function Radix2FFT{T}(fftsize::Integer) where T<:AbstractFloat # type-stability ✓
        # fftsize should be power of 2
        cache   = Vector{Complex{T}}(undef, fftsize)
        twiddle = Vector{Complex{T}}(undef, fftsize >> 1)

        return new{T}(cache, twiddle!(twiddle), fftsize, isone(pwr2(fftsize) & 1))
    end

    Radix2FFT(fftsize::Integer) = Radix2FFT{Float64}(fftsize)
end

"""
    fft!(x::AbstractVector{Complex{T}}, f::Radix2FFT{T}) where T<:AbstractFloat

Perform Fast-Fourier Transform (FFT) in place on the time-domain signal sequence `x` using an FFT kernel `f`.
The input sequence `x` will be overwritten by the FFT result.
"""
function fft!(x::VecI{Complex{T}}, f::Radix2FFT{T}) where T<:AbstractFloat
    cache = f.cache
    if f.ifswap
        @simd for i in eachindex(cache)
            @inbounds cache[i] = x[i]
        end
        ditnn!(cache, x, f.twiddle, f.fftsize >> 1)
    else
        ditnn!(x, cache, f.twiddle, f.fftsize >> 1)
    end

    return x
end

"""
    fft(x::AbstractVector{Complex{T}}, f::Radix2FFT{T}) where T<:AbstractFloat
    fft(x::AbstractVector{T},          f::Radix2FFT{T}) where T<:AbstractFloat

Similar to `fft!`, the function will pass the **copy** of the given signal sequence to the `fft!` routine.
"""
fft(x::VecI{Complex{T}}, f::Radix2FFT{T}) where T<:AbstractFloat = fft!(copy(x), f)

function fft(x::VecI{T}, f::Radix2FFT{T}) where T<:AbstractFloat
    cx = similar(x, Complex{T})
    @simd for i in eachindex(cx)
        @inbounds cx[i] = x[i]
    end
    return fft!(cx, f)
end

"""
    ifft!(x::AbstractVector{Complex{T}}, f::Radix2FFT{T}) where T<:AbstractFloat

Perform Inverse-Fast-Fourier Transform (IFFT) in place on the frequency-domain signal sequence `x` using an FFT kernel `f`.
The input sequence `x` will be overwritten by the IFFT result.
"""
function ifft!(x::VecI{Complex{T}}, f::Radix2FFT{T}) where T<:AbstractFloat
    fftsize = f.fftsize
    cache = f.cache

    if f.ifswap
        @simd for i in eachindex(x)
            @inbounds cache[i] = conj(x[i])
        end
        ditnn!(cache, x, f.twiddle, fftsize >> 1)
    else
        @simd for i in eachindex(x)
            @inbounds x[i] = conj(x[i])
        end
        ditnn!(x, cache, f.twiddle, fftsize >> 1)
    end

    @simd for i in eachindex(x)
        @inbounds x[i] = conj(x[i]) / fftsize
    end

    return x
end

"""
    ifft(x::AbstractVector{Complex{T}}, f::Radix2FFT{T}) where T<:AbstractFloat
    ifft(x::AbstractVector{T},          f::Radix2FFT{T}) where T<:AbstractFloat

Similar to `ifft!`, the function will pass the **copy** of the given signal sequence to the `ifft!` routine.
"""
ifft(x::VecI{Complex{T}}, f::Radix2FFT{T}) where T<:AbstractFloat = ifft!(copy(x), f)

function ifft(x::VecI{T}, f::Radix2FFT{T}) where T<:AbstractFloat
    cx = similar(x, Complex{T})
    @simd for i in eachindex(cx)
        @inbounds cx[i] = x[i]
    end
    return fft!(cx, f)
end
