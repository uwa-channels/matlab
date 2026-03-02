function unpacked_channel = unpack(fs, array_index, channel, varargin)
% UNPACK Unpack the channel coefficients.
%
% UNPACKED_CHANNEL = UNPACK(FS, ARRAY_INDEX, CHANNEL) returns the
% unpacked channel impulse response resampled to FS along the time axis.
%
% UNPACKED_CHANNEL = UNPACK(FS, ARRAY_INDEX, CHANNEL, BUF_LEFT, BUF_RIGHT)
% specifies the buffer ratios for padding the delay axis.
%
% The function handles three tracking modes:
%   theta_hat (phase tracking only):
%     h_hat already contains drifting taps. Only phase is re-inserted
%     via exp(+j*theta_hat). No delay interpolation is performed.
%   phi_hat (delay tracking):
%     h_hat is static. Both phase and delay drift are re-inserted.
%     Delay drift: drift = phi_hat / (2*pi*fc).
%   f_resamp:
%     Adds a linear phase ramp for bulk Doppler. Its delay contribution
%     is always re-inserted via interpolation.
%
% Inputs:
%    fs                  - Desired sampling rate along time axis.
%    array_index         - Indices of the hydrophone.
%    channel             - Struct containing parameters and impulse responses.
%    varargin            - Optional parameters for buffers (padding).
%
% Outputs:
%    unpacked_channel    - Resampled output signal.
%
% Example:
%    See example_unpack.m.
%
% Other m-files required: None
% Subfunctions: None
% Toolbox required: Signal Processing Toolbox (for resample function).
% MAT-files required: Channel MAT-file.
%
% See also: replay.m
%
% Author: Zhengnan Li
% Email : uwa-channels@ofdm.link
%
% License: MIT
%
% Revision history:
%   - Apr.  1, 2025: Initial release.
%   - Feb. 27, 2026: Fixed delay insertion for theta_hat vs phi_hat.
%

%% Unpacking variables
params = channel.params;
h_hat = channel.h_hat(:, array_index, :);

%% Parameters
fs_delay = params.fs_delay;
fs_time = params.fs_time;
fs_time_desired = fs;
fc = params.fc;
K = size(h_hat, 1);
M = size(h_hat, 2);
T = size(h_hat, 3);

N_phi = ceil(T * fs_delay / fs_time);
phase_all = zeros(length(array_index), N_phi);
phase_drift = zeros(length(array_index), N_phi);

if isfield(channel, 'theta_hat')
  phase_all = phase_all + channel.theta_hat(array_index, :);
end

if isfield(channel, 'phi_hat')
  phase_all = phase_all + channel.phi_hat(array_index, :);
  phase_drift = phase_drift + channel.phi_hat(array_index, :);
end

if isfield(channel, 'f_resamp')
  f_resamp_phase = (1 / channel.f_resamp - 1) * 2 * pi * fc ...
    * (1:N_phi) / fs_delay;
  phase_all = phase_all + f_resamp_phase;
  phase_drift = phase_drift + f_resamp_phase;
end

%% Allocate buffer
if nargin == 3
  buffer_left = 0.1;
  buffer_right = 0.1;
else
  buffer_left = varargin{1};
  buffer_right = varargin{2};
end

h_hat = [zeros(ceil(K*buffer_left), M, T); h_hat; zeros(ceil(K*buffer_right), M, T)];
K = size(h_hat, 1);

%% Sample rate conversion
[p1, q1] = rat(fs_time_desired / fs_time);
[p2, q2] = rat(fs_time_desired / fs_delay);
delays = (0:K-1) ./ fs_delay;

%% Unpack
has_phase = isfield(channel, 'theta_hat') || isfield(channel, 'phi_hat') || isfield(channel, 'f_resamp');
has_drift = isfield(channel, 'phi_hat') || isfield(channel, 'f_resamp');

unpacked_channel = zeros(K, length(array_index), ceil(T*p1/q1));

for m = 1:length(array_index)
  h_hat_m = squeeze(h_hat(:, m, :));
  h_resampled = resample(h_hat_m, p1, q1, 'Dimension', 2);

  if has_phase
    phase_resampled = resample(phase_all(m, :), p2, q2);
    h_resampled = h_resampled .* exp(1j * phase_resampled);
  end

  if has_drift
    drift_resampled = resample(phase_drift(m, :), p2, q2);
    drift = drift_resampled ./ (2 * pi * fc);
    for t = 1:size(h_resampled, 2)
      h_resampled(:, t) = interp1(delays, h_resampled(:, t), ...
        delays + drift(t), 'spline');
    end
  end

  unpacked_channel(:, m, :) = h_resampled;
end

unpacked_channel = unpacked_channel ./ max(abs(unpacked_channel), [], 'all');

end

% [EOF]
