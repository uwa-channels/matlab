# UWA-Channels MAT File Format Specification

The `uwa-channels`-compatible `.mat` files **must** be saved with the following flags:

* `-v7.3` to support large variables and `HDF5` file format.
* `-nocompression` to accelerate loading speed.

## Required Fields

Each `.mat` file must include the following variables:

### `h_hat`

* **Type**: Multi-dimensional numeric tensor.
* **Dimensions**: `[delay, receiver, time]`
* **Units**:

  * Delay axis: seconds
  * Time axis: seconds
  * Amplitude: complex baseband impulse response (unitless)
* **Description**: The estimated time-varying channel impulse response (TVIR).

### `params`

A MATLAB structure with the following scalar fields:

* `fs_delay` (Hz): Sampling rate along the *delay* axis. Determines the time resolution of the impulse response.
* `fs_time` (Hz): Sampling rate along the *time* axis. Determines how often the impulse response is sampled.
* `fc` (Hz): Center frequency of the signal used during channel estimation.

### `version`

* **Type**: Numeric scalar.
* **Description**: Dataset format version number.

## Optional Fields

### `theta_hat`

* **Type**: Numeric matrix of size `[receiver, time]`
* **Units**: Radians (interpreted as phase rotation angles)
* **Sampling Rate**: Must match `params.fs_delay`
* **Duration Constraint**: The third dimension of `h_hat` and the time dimension of `theta_hat` must span the same time duration, i.e.,
  ```
  size(theta_hat, 2) / params.fs_delay == size(h_hat, 3) / params.fs_time
  ``` 
* **Description**: For each hydrophone (receiver), `theta_hat` represents a time-varying *resampling factor* to be applied to the signal during or after convolution. The baseband received signal $v(t) is modeled as:

  $$
  v(t) = \sum_n d(n) h(t - nT) e^{j \theta(t)} + z(t)
  $$

  where $d(n)$ is the data symbol, $h(t)$ is the time-varying impulse response, $T$ is the symbol duration, and $\theta(t)$ is the phase correction term at time $t$.

### `f_resamp`

* **Type**: Scalar (double precision)
* **Units**: Unitless resampling factor
* **Description**: A *time-invariant* resampling factor applied to the output signal. This is typically the *inverse* of a resampling operation applied to remove the nominal Doppler frequency offset before channel estimation. It will be applied after the time-varying convolution.

### `meta`

The `meta` structure is optional but encouraged for documenting dataset context. The following fields are supported:

* `meta.experiment_name` (string): Name of the experiment.
* `meta.experiment_info` (string or struct): Descriptive details or notes about the experiment setup.
* `meta.array_info` (string or struct): Information about the hydrophone array (geometry, number of elements, etc.).
* `meta.authorship` (string or struct): Information about the data contributors and affiliations.

Users are free to add *any number of additional fields* to `meta` as needed to capture experiment-specific metadata.

## Dataset Configurations

The following configurations are currently supported:

* `h_hat` only: Supported, but not used in official datasets.
* `h_hat` and `theta_hat`: Used in *blue*, *red*, *green*, *purple*, *yellow*, *abyssal*, and *namikaze* datasets.
* `h_hat` and `f_resamp`: Used in *Watermark* datasets.
* `h_hat`, `theta_hat`, and `f_resamp`: Used in the *Sky* dataset.


