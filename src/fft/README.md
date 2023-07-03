# Fast-Fourier Transform (FFT) Module

This module provides a fast-Fourier transform (FFT) module for efficient FFT computations on time-domain and frequency-domain signal sequences. Currently, the module supports radix-2 FFT with a specific `T<:AbstractFloat` type. In the near future, it will also offer complex-to-complex FFT with arbitrary sequence sizes and more efficient real-to-complex FFT than complex-to-complex FFT for both radix-2 and arbitrary sequence sizes.

## # Instal and Import Module

To use this module, you can import this module by running the following command in the `Julia` REPL:

```julia
using Pkg
Pkg.add(url="https://github.com/FemtoPhysics/SpectroDSP.jl.git")

import SpectroDSP: FFT
```

## # Usage

For example, if you have a signal sequence `signal` of size `32`.

### ## Prepare Radix-2 Kernel

To initialize and obtain a kernel for radix-2 FFT, you can use the following constructor:

```julia
fft32 = FFT.Radix2FFT{eltype(signal)}(length(signal))
```

Or the default constructor with `Float64` type:

```julia
fft32 = FFT.Radix2FFT(length(signal))
```

### ## In-place and Copy-based FFT

```julia
# In-place FFT
spectrum = similar(signal, Complex{eltype(signal)})
FFT.fft!(spectrum, fft32)

# or

# Copy-based FFT
spectrum = FFT.fft(signal, fft32)
```

### ## In-place and Copy-based Inverse FFT (IFFT)

```julia
# In-place IFFT
recovered_signal = similar(spectrum, Complex{eltype(signal)})
FFT.ifft!(recovered_signal, fft32)

# or

# Copy-based IFFT
FFT.recovered_signal = ifft(spectrum, fft32)
```

## # Future Features

Soon, this module will provide additional functionality:

1. Complex-to-complex FFT with arbitrary sequence sizes.
2. More efficient real-to-complex FFT (compared to complex-to-complex FFT) for both radix-2 and arbitrary sequence sizes.

### TODO

- [x] C2C-FFT for radix-2 sizes
- [ ] R2C-FFT for radix-2 sizes
- [ ] C2C-FFT for arbitrary sizes
- [ ] R2C-FFT for arbitrary sizes

Stay tuned for updates and new releases!