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
%     is always re-inserted via interpolation. A scalar applies to every
%     hydrophone; a vector holds one factor per hydrophone in the file and
%     is indexed by ARRAY_INDEX.
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
% Subfunctions: validate_inputs
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
%   - Apr. 20, 2025: Added support for channels without phase tracking;
%                    accumulated phase contributions from theta_hat and
%                    f_resamp separately for Watermark compatibility.
%   - Feb. 27, 2026: Separated phi_hat (delay + phase) from theta_hat
%                    (phase only); restructured phase and drift fields.
%   - Jun. 10, 2026: Removed global normalization; fixed spline
%                    extrapolation by clamping out-of-range values to zero;
%                    replaced resample with interp1 (pchip) for phase and
%                    drift to eliminate FIR edge-effect ripples.
%   - Sep. 16, 2026: Fixed the origin of the f_resamp phase ramp.  It was
%                    built on (1:N_phi) while t_orig, the grid it is then
%                    interpolated from, is (0:N_phi-1)/fs_delay, so the ramp
%                    started one sample of fs_delay in.  That put a constant
%                    phase rotation, and through phase_drift a constant delay
%                    offset, on every unpacked tap.
%   - Oct. 27, 2026: phi_hat takes precedence over theta_hat when a file
%                    carries both, as in replay.m; their phases used to add.
%   - Nov.  3, 2026: f_resamp may be a vector with one factor per
%                    hydrophone in the file, indexed by array_index; each
%                    hydrophone gets its own phase ramp and delay drift.
%                    Added validate_inputs, which checks the channel
%                    version, fs, array_index, the buffers and f_resamp
%                    before unpacking, as replay.m does.
%

%% Buffers and simple checks
if nargin == 3
    buffer_left = 0.1;
    buffer_right = 0.1;
else
    assert(nargin == 5, 'unpack:invalidBuffer', ...
        'Give both buffers, BUF_LEFT and BUF_RIGHT, or neither.');
    buffer_left = varargin{1};
    buffer_right = varargin{2};
end
validate_inputs(fs, array_index, channel, buffer_left, buffer_right);

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

% phi_hat is the whole tracked phase (est.m saves its one trajectory as
% phi_hat with delay tracking and as theta_hat without), so it takes
% precedence over theta_hat, as in replay.m.
if isfield(channel, 'phi_hat')
    phase_all = phase_all + channel.phi_hat(array_index, :);
    phase_drift = phase_drift + channel.phi_hat(array_index, :);
elseif isfield(channel, 'theta_hat')
    phase_all = phase_all + channel.theta_hat(array_index, :);
end

if isfield(channel, 'f_resamp')
    % One factor per hydrophone in the file, or one for all of them.
    f_resamp = channel.f_resamp(:);
    if ~isscalar(f_resamp)
        f_resamp = f_resamp(array_index);
    end
    % Sample n of the trajectory is at (n-1)/fs_delay, the grid t_orig below
    % is built on, so the ramp has to start at zero and not at one sample in.
    f_resamp_phase = (1 ./ f_resamp - 1) * 2 * pi * fc ...
        * (0:N_phi-1) / fs_delay;
    phase_all = phase_all + f_resamp_phase;
    phase_drift = phase_drift + f_resamp_phase;
end

%% Allocate buffer
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

t_orig = (0:N_phi-1) ./ fs_delay;

for m = 1:length(array_index)
    h_hat_m = squeeze(h_hat(:, m, :));
    h_resampled = resample(h_hat_m, p1, q1, 'Dimension', 2);
    t_new = (0:size(h_resampled, 2)-1) ./ fs_time_desired;

    if has_phase
        phase_resampled = interp1(t_orig, phase_all(m, :), t_new, 'pchip', 'extrap');
        h_resampled = h_resampled .* exp(1j * phase_resampled);
    end

    if has_drift
        drift_resampled = interp1(t_orig, phase_drift(m, :), t_new, 'pchip', 'extrap');
        drift = drift_resampled ./ (2 * pi * fc);
        for t = 1:size(h_resampled, 2)
            h_resampled(:, t) = interp1(delays, h_resampled(:, t), ...
                delays + drift(t), 'spline', 0);
        end
    end

    unpacked_channel(:, m, :) = h_resampled;
end

end


function validate_inputs(fs, index, channel, buffer_left, buffer_right)
% Channel version checking
assert(channel.version >= 1.0, 'unpack:invalidVersion', ...
    'The minimum version of the channel matrix is 1.0, and you have %.1f.', ...
    channel.version);

% Check the sampling rate
assert(isnumeric(fs) && isscalar(fs) && isreal(fs) && isfinite(fs) && fs > 0, ...
    'unpack:invalidFs', 'fs must be a positive, finite scalar.');

% Check the hydrophone indices
N = size(channel.h_hat, 2);
assert(~isempty(index) && all(index == round(index)), ...
    'unpack:invalidArrayIndex', 'array_index must be positive integers.');
assert(length(unique(index)) == length(index), 'unpack:invalidArrayIndex', ...
    'array_index must not contain duplicate indices.');
assert(min(index) >= 1 && max(index) <= N, 'unpack:invalidArrayIndex', ...
    'array_index must be positive integers and cannot exceed %d.', N);

% Check the buffers
for b = {buffer_left, buffer_right}
    assert(isnumeric(b{1}) && isscalar(b{1}) && isreal(b{1}) ...
        && isfinite(b{1}) && b{1} >= 0, 'unpack:invalidBuffer', ...
        'The buffers must be non-negative, finite scalars.');
end

% Check the resampling factors
if isfield(channel, 'f_resamp')
    f = channel.f_resamp;
    assert(isnumeric(f) && isreal(f) && all(isfinite(f(:))) && all(f(:) > 0), ...
        'unpack:invalidFResamp', 'f_resamp must be positive and finite.');
    assert(numel(f) == 1 || numel(f) == N, 'unpack:invalidFResamp', ...
        'f_resamp must be a scalar or have one factor per hydrophone (%d), not %d.', ...
        N, numel(f));
end
end

% [EOF]
