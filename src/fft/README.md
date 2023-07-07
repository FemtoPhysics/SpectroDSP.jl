# Fast-Fourier Transform (FFT) Module

This module provides a fast-Fourier transform (FFT) module for efficient FFT computations on time-domain and frequency-domain signal sequences. Currently, the module supports radix-2 FFT and arbitrary sequence sizes with a specific `T<:AbstractFloat` type. In the near future, it will also offer more efficient real-to-complex FFT than complex-to-complex FFT for both radix-2 and arbitrary sequence sizes.

## # Instal and Import Module

To use this module, you can import this module by running the following command in the `Julia` REPL:

```julia
using Pkg
Pkg.add(url="https://github.com/FemtoPhysics/SpectroDSP.jl.git")

import SpectroDSP: FFT
```

## # Usage

### ## Prepare FFT Kernel

#### Prepare Radix-2 Kernel

For example, if you have a signal sequence `signal` of size `32`.

To initialize and obtain a kernel for radix-2 FFT, you can use the following constructor:

```julia
fftkernel = FFT.Radix2FFT{eltype(signal)}(length(signal))
```

Or the default constructor with `Float64` type:

```julia
fftkernel = FFT.Radix2FFT(length(signal))
```

#### Prepare Non-Radix-2 Kernel

For example, if you have a signal sequence `signal` of size `28`.

To initialize and obtain a kernel for radix-2 FFT, you can use the following constructor:

```julia
fftkernel = FFT.NonRadix2FFT{eltype(signal)}(length(signal))
```

Or the default constructor with `Float64` type:

```julia
fftkernel = FFT.NonRadix2FFT(length(signal))
```

### ## In-place and Copy-based FFT

```julia
# In-place FFT
spectrum = similar(signal, Complex{eltype(signal)})
FFT.fft!(spectrum, fftkernel)

# or

# Copy-based FFT
spectrum = FFT.fft(signal, fftkernel)
```

### ## In-place and Copy-based Inverse FFT (IFFT)

```julia
# In-place IFFT
recovered_signal = similar(spectrum, Complex{eltype(signal)})
FFT.ifft!(recovered_signal, fftkernel)

# or

# Copy-based IFFT
recovered_signal = FFT.ifft(spectrum, fftkernel)
```

### ## Obtain FFT Frequencies

```julia
freq = FFT.fftfreq(length(signal), Δt) # Δt is the timestep
```

### ## Obtain FFT Amplitude Spectrum

```julia
ampl = FFT.fftampl!(similar(signal), spectrum)
```

## # Future Features

Soon, this module will provide additional functionality:

1. Complex-to-complex FFT with arbitrary sequence sizes.
2. More efficient real-to-complex FFT (compared to complex-to-complex FFT) for both radix-2 and arbitrary sequence sizes.

### TODO

- [x] C2C-FFT for radix-2 sizes
- [ ] R2C-FFT for radix-2 sizes
- [x] C2C-FFT for arbitrary sizes
- [ ] R2C-FFT for arbitrary sizes
- [x] C2C-IFFT for radix-2 sizes
- [ ] C2C-IFFT for arbitrary sizes
- [x] FFT Frequency
- [x] FFT Amplitude

Stay tuned for updates and new releases!