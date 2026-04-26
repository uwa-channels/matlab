[![CI](https://github.com/uwa-channels/matlab/actions/workflows/ci.yml/badge.svg)](https://github.com/uwa-channels/matlab/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/uwa-channels/matlab/graph/badge.svg?token=NQ1M28NGYM)](https://codecov.io/gh/uwa-channels/matlab)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19643731.svg)](https://doi.org/10.5281/zenodo.19643731)


# Underwater Acoustic Channel Toolbox — MATLAB / Octave

[![Generic badge](https://img.shields.io/badge/MATLAB-R2021a-BLUE.svg)](https://shields.io/) with [Signal Processing Toolbox™](https://www.mathworks.com/products/signal.html) or [![Generic badge](https://img.shields.io/badge/Octave-9.0-BLUE.svg)](https://shields.io) with the [signal](https://gnu-octave.github.io/packages/signal/) and [statistics](https://gnu-octave.github.io/packages/statistics/) packages.

MATLAB®/Octave toolbox for replaying signals through measured underwater acoustic channels, generating realistic ocean noise, and unpacking stored impulse responses.  To learn more about the channels, check out the [documentation](https://uwa-channels.github.io/).

Please report bugs and suggest enhancements by [creating a new issue](https://github.com/uwa-channels/matlab/issues).  We welcome your feedback.  See [CONTRIBUTING.md](CONTRIBUTING.md) for more information.

## Functions

| Function | Description |
|----------|-------------|
| `replay` | Pass a passband signal through a measured underwater acoustic channel. |
| `noisegen` | Generate realistic ocean noise: pink Gaussian (17 dB/decade), spatially-correlated Gaussian, or impulsive (symmetric α-stable). |
| `unpack` | Reconstruct the full time-varying impulse response from the compressed representation. |

## Quick start

Download the channel MAT-files from [here](https://zenodo.org/records/19643731) and place them where MATLAB®/Octave can find them.

### Replay and noise generation

```matlab
channel = load('blue_1.mat');
noise = load('blue_1_noise.mat');

y = replay(input, fs, array_index, channel);
w = noisegen(size(y), fs, array_index, noise);
r = y + 0.05 * w;
```

See `examples/example_replay.m` for a complete example that generates a BPSK signal, replays it through the `blue_1` channel, adds noise, and plots the received signal, cross-correlation, and spectrum.

### Unpack

```matlab
channel = load('blue_1.mat');
unpacked = unpack(fs_time, array_index, channel);
```

See `examples/example_unpack.m` for details.

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=uwa-channels/matlab)

## Channel format

Each channel MAT-file contains:

| Variable | Description |
|----------|-------------|
| `h_hat` | Estimated impulse response, size (K × M × T) |
| `theta_hat` or `phi_hat` | Phase or delay-phase trajectory, size (M × N) |
| `params` | Struct with `fs_delay`, `fs_time`, `fc` |
| `meta` | Estimation metadata (see [estimate](https://github.com/uwa-channels/estimate) repo) |
| `version` | File format version |

Each noise MAT-file contains:

| Field | Description |
|-------|-------------|
| `Fs` | Sampling rate at which noise statistics were measured [Hz] |
| `R` | Signal bandwidth [Hz] |
| `alpha` | Stability index (2 = Gaussian, < 2 = impulsive) |
| `beta` | Mixing coefficients, size (M × M × K) |
| `fc` | Center frequency [Hz] |
| `rms_power` | Per-channel RMS power scaling, size (M × 1) |
| `version` | Noise struct version |

## Tests

This repository includes automated testing via [GitHub Actions](https://github.com/uwa-channels/matlab/actions).  The [tests](/tests) folder contains three test suites:

| Test | What it verifies |
|------|-----------------|
| `testReplay` | Generates random mobile channels ({static, mobile} × {theta\_hat, phi\_hat}), transmits a signal, and checks that cross-correlation peaks match the true multipath structure. |
| `testNoise` | Verifies output size, spectral shape (17 dB/decade), spatial correlation (theoretical vs. sample), bandpass filtering, rms\_power scaling, Gaussianity (α = 2), and heavy-tail behavior (α < 2). |
| `testUnpack` | Tests all tracking modes (none, theta\_hat, phi\_hat, f\_resamp, and combinations) for correct impulse response reconstruction. |

Tests run automatically on every push, ensuring continued correctness of the core functions.

## Related repositories

| Repository | Description |
|------------|-------------|
| [uwa-channels/estimate](https://github.com/uwa-channels/estimate) | Channel estimation from single-carrier signals, with visualization. |
| [uwa-channels/python](https://github.com/uwa-channels/python) | Python implementation of the replay toolbox. |

## License

The license is available in the [LICENSE](LICENSE) file within this repository.

© 2025–2026, Underwater Acoustic Channels Group.
