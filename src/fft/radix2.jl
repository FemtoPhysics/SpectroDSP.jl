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

    Radix2FFT{T}(fftsize::Int64) where T<:AbstractFloat

Initialize and return a kernel with `fftsize` and default type `Float64`.

    Radix2FFT(fftsize::Int64)

The value of `fftsize` is strictly restricted to be 2 to an integer power.
"""
struct Radix2FFT{T<:AbstractFloat} <: FFTKernel{T}
    cache  ::Vector{Complex{T}}
    twiddle::Vector{Complex{T}}
    fftsize::Int
    ifswap ::Bool

    function Radix2FFT{T}(fftsize::Int64) where T<:AbstractFloat # type-stability ✓
        iszero(clp2(fftsize) - fftsize) || throw(DomainError(fftsize))
        cache   = Vector{Complex{T}}(undef, fftsize)
        twiddle = Vector{Complex{T}}(undef, fftsize >> 1)

        return new{T}(cache, twiddle!(twiddle), fftsize, isone(pwr2(fftsize) & 1))
    end

    Radix2FFT(fftsize::Int64) = Radix2FFT{Float64}(fftsize)
end

"""
Non-Radix2 Kernel for 1D FFT

Initialize and return a kernel with `fftsize` and specific type `T<:AbstractFloat`.

    NonRadix2FFT{T}(fftsize::Int64) where T<:AbstractFloat

Initialize and return a kernel with `fftsize` and default type `Float64`.

    NonRadix2FFT(fftsize::Int64)

The value of `fftsize` is strictly restricted **NOT** to be 2 to an integer power.
"""
struct NonRadix2FFT{T<:AbstractFloat} <: FFTKernel{T}
    cache0   ::Vector{Complex{T}} # Cache for Radix-2 ops
    cache1   ::Vector{Complex{T}} # Cache for h, diag(h)
    cache2   ::Vector{Complex{T}} # Cache for y, Y, Z, z
    twiddle  ::Vector{Complex{T}} # Twddles for Radix-2 ops
    circulant::Vector{Complex{T}} # Chirps for h
    fftsize  ::Int
    extsize  ::Int
    ifswap   ::Bool

    function NonRadix2FFT{T}(fftsize::Int64) where T<:AbstractFloat # type-stability ✓
        iszero(clp2(fftsize) - fftsize) && throw(DomainError(fftsize))
        extsize   = clp2((fftsize - 1) << 1) # extended size
        cache0    = Vector{Complex{T}}(undef, extsize)
        cache1    = Vector{Complex{T}}(undef, extsize)
        cache2    = Vector{Complex{T}}(undef, extsize)
        twiddle   = Vector{Complex{T}}(undef, extsize >> 1)
        circulant = Vector{Complex{T}}(undef, extsize)
        ifswap    = isone(pwr2(extsize) & 1)

        return new{T}(
            cache0, cache1, cache2, twiddle!(twiddle),
            circulantChirp!(circulant, extsize, fftsize),
            fftsize, extsize, isone(pwr2(extsize) & 1)
        )
    end

    NonRadix2FFT(fftsize::Int64) = NonRadix2FFT{Float64}(fftsize)
end

# = = = = = = = = = = = = = = = = = = = = = #
# Foward Radix-2 Fast-Fourier Transform     #
# = = = = = = = = = = = = = = = = = = = = = #

function fft!(signal::VecI{Complex{T}}, buffer::VecI{Complex{T}}, twiddle::VecI{Complex{T}}, fftsize::Int, ifswap::Bool) where T<:AbstractFloat
    if ifswap
        @simd for i in eachindex(buffer)
            @inbounds buffer[i] = signal[i]
        end
        ditnn!(buffer, signal, twiddle, fftsize >> 1)
    else
        ditnn!(signal, buffer, twiddle, fftsize >> 1)
    end

    return signal
end

# = = = = = = = = = = = = = = = = = = = = = #
# Backward Radix-2 Fast-Fourier Transform   #
# = = = = = = = = = = = = = = = = = = = = = #

function ifft!(signal::VecI{Complex{T}}, buffer::VecI{Complex{T}}, twiddle::VecI{Complex{T}}, fftsize::Int, ifswap::Bool) where T<:AbstractFloat
    if ifswap
        @simd for i in eachindex(signal)
            @inbounds buffer[i] = conj(signal[i])
        end
        ditnn!(buffer, signal, twiddle, fftsize >> 1)
    else
        @simd for i in eachindex(signal)
            @inbounds signal[i] = conj(signal[i])
        end
        ditnn!(signal, buffer, twiddle, fftsize >> 1)
    end

    @simd for i in eachindex(signal)
        @inbounds signal[i] = conj(signal[i]) / fftsize
    end

    return signal
end

"""
    fft!(x::AbstractVector{Complex{T}}, f::FFTKernel{T}) where T<:AbstractFloat

Perform Fast-Fourier Transform (FFT) in place on the time-domain signal sequence `x` using an FFT kernel `f`.
The input sequence `x` will be overwritten by the FFT result.
"""
fft!(x::VecI{Complex{T}}, f::Radix2FFT{T}) where T<:AbstractFloat = fft!(x, f.cache, f.twiddle, f.fftsize, f.ifswap)

function fft!(x::VecI{Complex{T}}, f::NonRadix2FFT{T}) where T<:AbstractFloat
    fftsize   = f.fftsize
    extsize   = f.extsize
    cache0    = f.cache0
    cache1    = f.cache1
    cache2    = f.cache2
    twiddle   = f.twiddle
    circulant = f.circulant
    ifswap    = f.ifswap

    copyto!(cache1, circulant)

    # diag(h) = FFT(h)
    fft!(cache1, cache0, twiddle, extsize, ifswap)

    # extend signal to y
    @inbounds for i in 1:fftsize
        cache2[i] = x[i] / circulant[i]
    end
    @inbounds for i in fftsize+1:extsize
        cache2[i] = complex(0.0, 0.0)
    end

    # Y = FFT(y)
    fft!(cache2, cache0, twiddle, extsize, ifswap)

    # Scaling: Z = h .* Y
    @inbounds for i in eachindex(cache2)
        cache2[i] *= cache1[i]
    end

    # z = IFFT(Z)
    ifft!(cache2, cache0, twiddle, extsize, ifswap)

    # Reconstruction
    @inbounds for i in eachindex(x)
        x[i] = cache2[i] / circulant[i]
    end

    return x
end

"""
    fft(x::AbstractVector{Complex{T}}, f::FFTKernel{T}) where T<:AbstractFloat
    fft(x::AbstractVector{T},          f::FFTKernel{T}) where T<:AbstractFloat

Similar to `fft!`, the function will pass the **copy** of the given signal sequence to the `fft!` routine.
"""
fft(x::VecI{Complex{T}}, f::FFTKernel{T}) where T<:AbstractFloat = fft!(copy(x), f)

function fft(x::VecI{T}, f::FFTKernel{T}) where T<:AbstractFloat
    cx = similar(x, Complex{T})
    @simd for i in eachindex(cx)
        @inbounds cx[i] = x[i]
    end
    return fft!(cx, f)
end

"""
    ifft!(x::AbstractVector{Complex{T}}, f::FFTKernel{T}) where T<:AbstractFloat

Perform Inverse-Fast-Fourier Transform (IFFT) in place on the frequency-domain signal sequence `x` using an FFT kernel `f`.
The input sequence `x` will be overwritten by the IFFT result.
"""
ifft!(x::VecI{Complex{T}}, f::Radix2FFT{T}) where T<:AbstractFloat = ifft!(x, f.cache, f.twiddle, f.fftsize, f.ifswap)

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
