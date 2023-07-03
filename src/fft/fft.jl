module FFT

const VecI = AbstractVector

## power of radix 2 (flooring mode)
function pwr2(x::Int64)
    x > 0 || error("pwr2(x): x should be positive!")
    r = 0
    x > 4294967295 && (x >>= 32; r += 32)
    x >      65535 && (x >>= 16; r += 16)
    x >        255 && (x >>=  8; r +=  8)
    x >         15 && (x >>=  4; r +=  4)
    x >          3 && (x >>=  2; r +=  2)
    x >          1 && (x >>=  1; r +=  1)

    return r
end

## power of radix 2 (ceiling mode)
function clp2(x::Int64)
    x ≡ 0 && return 1
    x ≡ 1 && return 2
    x = x - 1
    x = x | (x >> 1)
    x = x | (x >> 2)
    x = x | (x >> 4)
    x = x | (x >> 8)
    x = x | (x >> 16)
    return x + 1
end

include("./twiddle.jl")
include("./radix2.jl")

end
