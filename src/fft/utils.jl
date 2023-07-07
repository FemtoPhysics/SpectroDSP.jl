"""
    fftshift!(x::VecI)

Shift the zero-frequency component to the center of the spectrum.
"""
function fftshift!(x::VecI)
    N = length(x)
    M = N >> 1
    if iszero(N & 1)
        for i in 1:M
            swap!(x, i, i+m)
        end
    else
        i = j = 1
        s = @inbounds x[i]
        t = zero(eltype(x))
        for _ in 1:N
            j = i + M
            if j > N
                j -= N
            end
            @inbounds t = x[j]
            @inbounds x[j] = s
            s = t
            i = j
        end
    end
    return x
end

function fftfreq!(f::VecI, Δt::Real)
    N    = length(f)
    Δf   = inv(Δt * N)
    Np1  = N + 1
    Nby2 = N >> 1

    if iszero(N & 1)
        @inbounds for i in 1:Nby2
            f[i] = Δf * (i - 1)
        end
        @inbounds for i in Nby2+1:N
            f[i] = Δf * (i - Np1)
        end
    else
        @inbounds for i in 1:Nby2+1
            f[i] = Δf * (i - 1)
        end
        @inbounds for i in Nby2+2:N
            f[i] = Δf * (i - Np1)
        end
    end

    return f
end

"""
Discrete Fourier Transform sample frequencies

    fftfreq(n::Int, Δt::Real)

The returned array `f` contains the `n`-long frequencies with zero at the start.
`Δt` is the timestep of the equidistant time series.
"""
fftfreq(n::Int, Δt::Real) = fftfreq!(Vector{Float64}(undef, n), Δt)

"""
Discrete Fourier Transform amplitudes

    fftampl!(ampl::AbstractVector, spec::AbstractVector)

Compute the FFT real-valued amplitudes from a complex spectrum `spec` and store them in the vector `ampl`.
"""
function fftampl!(ampl::VecI, spec::VecI)
    hfsz = length(ampl) >> 1
    @inbounds for i in eachindex(ampl)
        ampl[i] = apy2(spec[i]) / hfsz
    end
    return ampl
end
