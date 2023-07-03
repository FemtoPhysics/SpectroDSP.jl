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
