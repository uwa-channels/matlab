# UWA-Channels MAT File Format Specification

This specification is also published at <https://uwa-channels.github.io/docs/file_specifications/>; keep the two copies in sync.

The `uwa-channels`-compatible `.mat` files **must** be saved with the following flags:

* `-v7.3` to support large variables and the `HDF5` file format.
* `-nocompression` to speed up loading.

The Python loaders `load_channel` and `load_noise` also read MAT v5 files, but files you distribute should be saved with `-v7.3 -nocompression`.

## Required fields

Each `.mat` file must include the following variables:

### `h_hat`

* **Type**: Multi-dimensional complex tensor.
* **Dimensions**: `[delay, receiver, time]`
* **Units**:
  * Delay axis: sampled at `params.fs_delay` [Hz], denoted $f_s = 1/T_s$
  * Time axis: sampled at `params.fs_time` [Hz]
  * Amplitude: complex baseband impulse response (unitless)
* **Description**: The estimated time-varying channel impulse response (TVIR) as a function of delay $\tau$ and time $t$. When delay tracking is used, this is the drift-free response $\hat{\underline{h}}(\tau, t)$; otherwise it is the drifting response $\hat{h}(\tau, t)$. For a given receiver and time index $n$, the slice along the delay axis is the channel vector $\hat{\underline{\mathbf{h}}}[n]$, or $\hat{\mathbf{h}}[n]$ in the drifting case.

### `params`

A structure with the following scalar fields:

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `fs_delay` | scalar | Hz | Sampling rate $f_s$ along the delay axis. |
| `fs_time` | scalar | Hz | Sampling rate along the time axis. |
| `fc` | scalar | Hz | Center frequency $f_c$ of the signal used during channel estimation. |

### `version`

* **Type**: Numeric scalar.
* **Description**: Dataset format version number (currently `1.0`). `replay` rejects files with a version below `1.0`.

## Phase/delay tracking fields (optional)

`phi_hat`, `theta_hat`, and `f_resamp` may coexist; `replay` and `unpack` apply them according to the following precedence:

* If none of the three is present, no Doppler correction is applied.
* `phi_hat` takes precedence over `theta_hat`: if `phi_hat` is present, it is used and `theta_hat` is ignored; if only `theta_hat` is present, it is used instead.
* `f_resamp` (see [Optional fields](#optional-fields)) is applied independently of the two phase/delay fields, whether or not either of them is present.

The estimation scripts save their one tracked phase as `phi_hat` with delay tracking and as `theta_hat` without, so a file normally carries only one of the two; `phi_hat` is the whole phase, which is why it wins.

### `phi_hat` (delay tracking)

* **Type**: Numeric matrix, size `[receiver, time]`
* **Units**: Radians
* **Sampling rate**: `params.fs_delay`
* **Description**: The phase estimate $\hat\varphi(nT_s)$, a time-varying phase that encodes *both* phase rotation and delay drift. In this mode, `h_hat` holds the drift-free response $\hat{\underline{h}}(\tau, t)$. Let $\Delta\tau$ denote the delay drift common to all propagation paths. It is related to the phase by

  $$\hat\varphi(nT_s) = -2\pi f_c\, \Delta\tau(nT_s)$$

  Unpacking first reinserts the phase, multiplying the drift-free response by $e^{j\hat\varphi(nT_s)}$ to obtain the signal $\bar y_{\text{r}}(nT_s)$. It then reinserts the delay drift by evaluating that signal at shifted instants,

  $$\bar v_{\text{r}}(nT_s) = \mathcal{I}\left[\bar y_{\text{r}}\left(nT_s + \frac{\hat\varphi(nT_s)}{2\pi f_c}\right)\right]$$

  where $\mathcal{I}[\cdot]$ denotes interpolation, implemented here as spline interpolation.

  `phi_hat` accumulates over the whole record, so its absolute value can reach hundreds of samples of delay. `replay` therefore reinserts it relative to the starting sample, $\hat\varphi(nT_s) - \hat\varphi(n_0T_s)$ for a replay that starts at sample $n_0$: the output starts aligned with the input wherever in the channel the replay starts, and it is padded by the largest delay the drift reaches within the replay.

### `theta_hat` (phase tracking only)

* **Type**: Numeric matrix, size `[receiver, time]`
* **Units**: Radians
* **Sampling rate**: `params.fs_delay`
* **Description**: The phase estimate $\hat\theta(nT_s)$, a time-varying phase correction. In this mode, `h_hat` contains the *drifting* impulse response $\hat{h}(\tau, t)$ (delay drift is embedded in the taps). Only the phase is tracked separately. The baseband received signal is modeled as:

  $$v(t) = \sum_n d(n)\, h(t - nT, t)\, e^{j\theta(t)} + w(t)$$

  where $d(n)$ are the transmitted data symbols, $h(\tau, t)$ is the time-varying impulse response at delay $\tau$ and time $t$ with drifting taps, $T$ is the symbol interval, $\theta(t)$ is the channel phase, and $w(t)$ is the additive complex baseband noise.

### Duration constraint

The time dimension of `theta_hat` or `phi_hat` and the third dimension of `h_hat` must span the same duration. The check is written below for `theta_hat`; the same applies to `phi_hat`:

```
size(theta_hat, 2) / params.fs_delay == size(h_hat, 3) / params.fs_time
```

## Optional fields

### `f_resamp`

* **Type**: Scalar (double precision)
* **Units**: Unitless resampling factor
* **Description**: A time-invariant resampling factor applied to the output signal. This is typically the inverse of a resampling operation applied to remove the nominal Doppler frequency offset before channel estimation. `replay` applies it after the time-varying convolution, by resampling the output. `unpack` instead converts it into a linear phase ramp $(1/f_{\text{resamp}} - 1)\,2\pi f_c t$, which it reinserts as both a phase and a delay drift, in the same way as `phi_hat`.

### `meta`

The `meta` structure is optional but strongly encouraged. The following fields are written by the [estimate](https://github.com/uwa-channels/estimate) scripts; `replay`, `noisegen`, and `unpack` do not read `meta`:

| Field | Type | Description |
|-------|------|-------------|
| `description` | string | Free-text description of the experiment. |
| `nsd` | scalar | Samples per symbol $N_s$ in the delay domain. |
| `nst` | scalar | Samples per symbol in the time domain. |
| `K_1` | scalar | Anti-causal filter length [symbols]. |
| `K_2` | scalar | Causal filter length [symbols]. |
| `fc` | scalar | Center frequency $f_c$ [Hz]. |
| `element_spacing` | scalar | Array element spacing $\ell$ [m]. |
| `vertical` | logical | `true` if vertical array. |
| `delay_tracking` | logical | `true` if delay tracking is enabled (`phi_hat` present). |
| `limit` | scalar | Lower dB limit for plotting. |
| `optim` | scalar | Optimizer used: `1` LMS, `2` RLS, `3` SFTF. |
| `mu` | scalar | LMS step size $\mu$ (when `optim == 1`). |
| `lambda` | scalar | Forgetting factor $\lambda$ (when `optim == 2` or `3`). |
| `regularization` | scalar | Regularization factor (when `optim == 2` or `3`). |
| `Kf_1` | scalar | PLL loop filter coefficient $K_{f_1}$. |
| `Kf_2` | scalar | PLL loop filter coefficient $K_{f_2}$. |
| `nslr` | scalar | Delay tracking rate (when `delay_tracking == true`). |
| `codename` | string | Short identifier for the channel (e.g., `"blue_1"`). |

Users are free to add additional fields to `meta` to capture experiment-specific metadata.

## Noise file format

Each noise `.mat` file contains the following fields:

| Field | Type | Description |
|-------|------|-------------|
| `Fs` | scalar | Sampling rate $f_s$ at which noise statistics were measured [Hz]. |
| `R` | scalar | Bandwidth $B$ of the recorded noise [Hz]. Despite the field name, this is not the symbol rate $R$ used on the channel pages. |
| `alpha` | scalar | Characteristic exponent $\alpha$ of the symmetric $\alpha$-stable distribution: `2` Gaussian, `< 2` impulsive. |
| `beta` | tensor `[M, M, L+1]` | Mixing coefficients $\beta_{ij}(kT_s)$ for spatiotemporal noise coloring, where $M$ is the number of array elements and $L$ is the maximum lag (both defined below). |
| `fc` | scalar | Center frequency $f_c$ [Hz]. |
| `version` | scalar | Noise struct version number. `noisegen` rejects files with a version below `1.0`. |

`noisegen` reads only `Fs`, `alpha`, `beta`, and `version`; it generates the noise at `Fs` and resamples it to the requested rate. `beta` is normalized over the whole array: the expected noise powers of all $M$ elements sum to $M$, and each element keeps its measured power relative to the others. The noise generated for a subset of elements (`array_index`) therefore need not sum to the size of the subset. When $\alpha < 2$ the variance is undefined, and the pseudo-power $2c^2$, with $c$ the scale parameter, takes its place. The noise generation function `noisegen` uses the mixing equation:

$$\hat{n}_i(nT_s) = \sum_{j=0}^{M-1}\sum_{k=0}^{L}\beta_{ij}(kT_s)\,\eta_j(nT_s + kT_s)$$

Here $\hat{n}_i$ is the synthetically generated noise on hydrophone $i$, with the hydrophone indices running as $i, j = 0, \ldots, M-1$. The sampling interval is $T_s = 1/f_s$, the sampling rate $f_s$ being stored in the `Fs` field, and $L$ is the maximum discrete time lag with non-negligible covariance, so that $\beta_{ij}(kT_s)$ has $L+1$ taps, $k = 0, \ldots, L$. The innovations $\eta_j(nT_s) \sim \mathcal{S}_\alpha(0, 1/\sqrt{2})$ are i.i.d. symmetric $\alpha$-stable variates with zero location and scale $1/\sqrt{2}$; this scale ensures that the distribution reduces to the standard Gaussian $\mathcal{N}(0, 1)$ when $\alpha = 2$. The sum runs forward in time, $\eta_j(nT_s + kT_s)$, which is the convention under which `beta` is estimated by [NoiseModels.jl](https://github.com/org-arl/NoiseModels.jl). For an array, mixing backward instead, $\eta_j(nT_s - kT_s)$, reverses the time direction of the correlation between elements while leaving the zero-lag covariance unchanged.
