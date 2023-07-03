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
