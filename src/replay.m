function output = replay(input, fs, array_index, channel, varargin)
% REPLAY Pass the signal through the given underwater acoustic channel.
%
% OUTPUT = REPLAY(INPUT, FS, ARRAY_INDEX, CHANNEL) returns the OUTPUT
% signal after being processed through the CHANNEL. FS is the sampling
% frequency of the INPUT signal, ARRAY_INDEX specifies the indices
% of the CHANNEL, and CHANNEL is a struct containing the estimated
% channel coefficients, digital phase-locked loop outputs, and
% other relevant parameters.
%
% OUTPUT = REPLAY(INPUT, FS, ARRAY_INDEX, CHANNEL, START) specifies
% the START index of the channel traces, for reproducibility.
%
% Inputs:
%    input               - Input signal (time-by-array).
%    fs                  - Sampling frequency of the input signal in Hz.
%    array_index         - Indices of the hydrophone.
%    channel             - Struct containing parameters and impulse responses.
%    varargin            - Optional parameter to specify start point in time.
%
% Outputs:
%    output              - Processed output signal.
%
% Example:
%    see example_replay.m.
%
% Other m-files required: None
% Subfunctions: validate_inputs, pwr
% Toolbox required: Signal Processing Toolbox (R) (resample function).
% MAT-files required: channel matfile.
%
% See also: noisegen.m
%
% Author: Zhengnan Li
% Email : uwa-channels@ofdm.link
%
% License: MIT
%
% Revision history:
%   - Apr. 1, 2025: initial release.
%
%

%% Simple checks
validate_inputs(input, fs, array_index, channel);

%% Unpacking variables
h_hat = channel.h_hat(:, array_index, :); % Channel matrix, delay x receiver x time.
fs_delay = channel.params.fs_delay; % Sampling rate in delay domain.
fs_time = channel.params.fs_time; % Sampling rate in time.
fc = channel.params.fc; % Channel center frequency.
M = length(array_index); % Number of array elements in user's input.
L = size(h_hat, 1); % Length of channel estimator.
if isfield(channel, 'theta_hat')
    theta_hat = channel.theta_hat(array_index, :); % PLL outputs, receiver x time.
end

%% Convert to baseband and resample the baseband to fs_delay
[p, q] = rat(fs_delay/fs);
baseband = input .* exp(-2j*pi*fc*(0:size(input, 1) - 1).'/fs);
baseband = resample(baseband, p, q);
T = length(baseband);

%% Assign random start point in time (for reproducibility only)
buffer = 20; % extra buffer for extrapolation
T_max = size(h_hat, 3) / fs_time * fs_delay;
if nargin ~= 5
    start = randi([0, T_max - T - L - buffer - 1])
else
    start = varargin{end}
end

%% Convolution and insert the drift
output = zeros(T+buffer+L, M);
baseband = [zeros(L-1, 1); baseband; zeros(L-1, 1)];
signal_time = ((0:T + L + buffer - 1) + start) ./ fs_delay;
signal_start = floor(min(signal_time)*fs_time);
signal_end = ceil(max(signal_time)*fs_time);
[p1, q1] = rat(fs_delay/fs_time);
for m = 1:M
    h_hat_m = flip(squeeze(h_hat(:, m, :)), 1);
    ir = resample(h_hat_m(:, signal_start:signal_end).', p1, q1, 'Dimension', 1);
    if isfield(channel, 'theta_hat')
        for t = 1:T + L - 1
            output(t, m) = ir(t, :) * baseband(t:t+L-1) .* exp(1j*theta_hat(m, t+start-1));
        end
        % Insert the drift
        drift = theta_hat(m, (0:T + L + buffer - 1)+start) ./ (2 * pi * fc);
        output(:, m) = interp1(signal_time, output(:, m), signal_time+drift, 'spline');
    else
        for t = 1:T + L - 1
            output(t, m) = ir(t, :) * baseband(t:t+L-1);
        end
    end
end

%% Resample to match the original sampling rate and upshift to fc
output_resampled = resample(output, q, p, 'Dimension', 1);
output_resampled = real(output_resampled.*exp(2j*pi*fc*(0:size(output_resampled, 1) - 1).'/fs));

if isfield(channel, 'f_resamp')
    [p2, q2] = rat(channel.f_resamp);
    output_resampled = resample(output_resampled, p2, q2, 'Dimension', 1);
end
output = 1 ./ sqrt(sum(pwr(output_resampled))) .* output_resampled;

end


function p = pwr(x)
p = mean(abs(x).^2, 1);
end


function validate_inputs(input, fs, index, channel)
% Channel version checking
assert(channel.version >= 1.0, ...
    'The minimum version of the channel matrix is 1.0, and you have %.1f.', ...
    channel.version);

% Get constants
[T, M] = size(input);
T = T / fs;
[~, N, T_max] = size(channel.h_hat);
T_max = T_max / channel.params.fs_time;

% Check the signal lengths
assert(T < T_max, ...
    'Duration of the input signal, %.2fms, should be no larger than %.2fms.', ...
    T*1e3, T_max*1e3);

% Check indices for two modes
assert(length(unique(index)) == length(index));
assert(max(index) <= N, ...
    'array_index must be positive integers and cannot exceed %d.', max(index));
assert(M <= N, 'The maximum supported number of channels is %d.', N);

end

% [EOF]
