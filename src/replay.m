function output = replay(input, fs, array_index, channel, varargin)
% REPLAY Passes a passband signal through an underwater acoustic channel.
%
% OUTPUT = REPLAY(INPUT, FS, ARRAY_INDEX, CHANNEL) returns the received OUTPUT
% signal after the INPUT passband signal is processed through the CHANNEL.
% FS is the sampling frequency of the INPUT signal, ARRAY_INDEX specifies the
% hydrophone indices, and CHANNEL is a struct containing the estimated
% channel coefficients and other relevant parameters.
%
% OUTPUT = REPLAY(INPUT, FS, ARRAY_INDEX, CHANNEL, START) specifies the START
% index of the channel traces for reproducibility.
%
% Inputs:
%    input               - Input passband signal (time-by-array).
%    fs                  - Sampling frequency of the input signal in Hz.
%    array_index         - Indices of the hydrophones.
%    channel             - Struct containing parameters and impulse responses.
%    varargin            - Optional parameter specifying the start time index.
%
% Outputs:
%    output              - Processed output signal.
%
% Example:
%    See example_replay.m.
%
% Other m-files required: None
% Subfunctions: validate_inputs
% Toolbox required: Signal Processing Toolbox (for resample function).
% MAT-files required: Channel MAT-file.
%
% See also: noisegen.m
%
% Author: Zhengnan Li
% Email : uwa-channels@ofdm.link
%
% License: MIT
%
% Revision history:
%   - Apr.  1, 2025: Initial release.
%   - May.  1, 2025: Replaced resample with interp1 (spline) for channel
%                    impulse response interpolation along the time axis.
%   - Mar.  1, 2026: Fixed start index off-by-one; distinguished phi_hat
%                    (delay + phase) from theta_hat (phase only).
%   - May. 30, 2026: Added power scaling to sqrt(M)/sqrt(sum(pwr)).
%   - Jul.  7, 2026: Removed normalization from replay.m.
%   - Jul.  8, 2026: Vectorized the time-varying convolution (sum over the
%                    L filter taps instead of a per-sample dot product).
%

%% Simple checks
validate_inputs(input, fs, array_index, channel);

%% Unpacking variables
fs_delay = channel.params.fs_delay; % Sampling rate in delay domain.
fs_time = channel.params.fs_time; % Sampling rate in time.
fc = channel.params.fc; % Channel center frequency.
M = length(array_index); % Number of array elements in user's input.
L = size(channel.h_hat, 1); % Length of channel estimator.

%% Convert to baseband and resample the baseband to fs_delay
[p, q] = rat(fs_delay/fs);
baseband = input .* exp(-2j*pi*fc*(0:size(input, 1) - 1).'/fs);
baseband = resample(baseband, p, q);
T = length(baseband);

%% Assign random start point in time (for reproducibility only)
buffer = 20; % extra buffer for extrapolation
T_max = size(channel.h_hat, 3) / fs_time * fs_delay;
if nargin ~= 5
    start = randi([1, T_max - T - L - buffer - 1]);
else
    start = varargin{end};
end

%% Convolution and insert the drift
output = zeros(T+buffer+L, M);
baseband = [zeros(L-1, 1); baseband; zeros(L-1, 1)];
channel_time = (0:size(channel.h_hat, 3) - 1) ./ fs_time;
signal_time = ((0:T + L + buffer - 1) + start) ./ fs_delay;
N = T + L - 1;
for m = 1:M
    h_hat_m = flip(squeeze(channel.h_hat(:, array_index(m), :)).', 2);
    ir = interp1(channel_time, h_hat_m, signal_time, 'spline', 0);

    % Time-varying convolution
    conv_out = zeros(N, 1);
    for l = 1:L
        conv_out = conv_out + ir(1:N, l) .* baseband(l:l+N-1);
    end

    if isfield(channel, 'phi_hat')
        phase = exp(1j * channel.phi_hat(array_index(m), start:start+N-1)).';
        output(1:N, m) = conv_out .* phase;
        % Insert the drift
        drift = channel.phi_hat(array_index(m), (0:T + L + buffer - 1)+start) ./ (2 * pi * fc);
        output(:, m) = interp1(signal_time, output(:, m), signal_time+drift, 'spline', 0);
    elseif isfield(channel, 'theta_hat')
        phase = exp(1j * channel.theta_hat(array_index(m), start:start+N-1)).';
        output(1:N, m) = conv_out .* phase;
    else
        output(1:N, m) = conv_out;
    end
end

%% Resample to match the original sampling rate and upshift to fc
output = resample(output, q, p, 'Dimension', 1);
output = 2 .* real(output .*exp(2j*pi*fc*(0:size(output, 1) - 1).'/fs));

if isfield(channel, 'f_resamp')
    [p2, q2] = rat(channel.f_resamp);
    output = resample(output, p2, q2, 'Dimension', 1);
end

end


function validate_inputs(input, fs, index, channel)
% Channel version checking
assert(channel.version >= 1.0, ...
    'The minimum version of the channel matrix is 1.0, and you have %.1f.', ...
    channel.version);

% Get constants
T = size(input, 1);
T = T / fs;
N = size(channel.h_hat, 2);
T_max = size(channel.h_hat, 3);
T_max = T_max / channel.params.fs_time;

% Check the signal lengths
assert(T < T_max, ...
    'Duration of the input signal, %.2fms, should be no larger than %.2fms.', ...
    T*1e3, T_max*1e3);

% Check the hydrophone indices
assert(length(unique(index)) == length(index), ...
    'array_index must not contain duplicate indices.');
assert(min(index) >= 1 && max(index) <= N, ...
    'array_index must be positive integers and cannot exceed %d.', N);
end

% [EOF]
