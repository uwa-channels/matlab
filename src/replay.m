function output = replay(input, fs, channel, varargin)
% REPLAY Passes a passband signal through an underwater acoustic channel.
%
% OUTPUT = REPLAY(INPUT, FS, CHANNEL) returns the received OUTPUT signal
% after the INPUT passband signal is processed through the CHANNEL.
% FS is the sampling frequency of the INPUT signal, and CHANNEL is a struct
% containing the estimated channel coefficients and other relevant parameters.
%
% OUTPUT = REPLAY(INPUT, FS, CHANNEL, START) specifies the START index of
% the channel traces for reproducibility.
%
% Inputs:
%    input               - Input passband signal (time-by-array).
%    fs                  - Sampling frequency of the input signal in Hz.
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
% Subfunctions: validate_inputs, pwr
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
%   - Jun. 10, 2026: Removed array_index parameter; function now uses all
%                    channels in channel.h_hat.
%

%% Simple checks
validate_inputs(input, fs, channel);

%% Unpacking variables
fs_delay = channel.params.fs_delay; % Sampling rate in delay domain.
fs_time = channel.params.fs_time; % Sampling rate in time.
fc = channel.params.fc; % Channel center frequency.
M = size(channel.h_hat, 2); % Number of array elements in user's input.
L = size(channel.h_hat, 1); % Length of channel estimator.

%% Convert to baseband and resample the baseband to fs_delay
[p, q] = rat(fs_delay/fs);
baseband = input .* exp(-2j*pi*fc*(0:size(input, 1) - 1).'/fs);
baseband = resample(baseband, p, q);
T = length(baseband);

%% Assign random start point in time (for reproducibility only)
buffer = 20; % extra buffer for extrapolation
T_max = size(channel.h_hat, 3) / fs_time * fs_delay;
if nargin ~= 4
    start = randi([1, T_max - T - L - buffer - 1])
else
    start = varargin{end}
end

%% Convolution and insert the drift
output = zeros(T+buffer+L, M);
baseband = [zeros(L-1, 1); baseband; zeros(L-1, 1)];
channel_time = (0:size(channel.h_hat, 3) - 1) ./ fs_time;
signal_time = ((0:T + L + buffer - 1) + start) ./ fs_delay;
for m = 1:M
    h_hat_m = flip(squeeze(channel.h_hat(:, m, :)).', 2);
    ir = interp1(channel_time, h_hat_m, signal_time, 'spline');
    if isfield(channel, 'phi_hat')
        for t = 1:T + L - 1
            output(t, m) = ir(t, :) * baseband(t:t+L-1) .* exp(1j*channel.phi_hat(m, t+start-1));
        end
        % Insert the drift
        drift = channel.phi_hat(m, (0:T + L + buffer - 1)+start) ./ (2 * pi * fc);
        output(:, m) = interp1(signal_time, output(:, m), signal_time+drift, 'spline');
    elseif isfield(channel, 'theta_hat')
        for t = 1:T + L - 1
            output(t, m) = ir(t, :) * baseband(t:t+L-1) .* exp(1j*channel.theta_hat(m, t+start-1));
        end
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
output = sqrt(M) ./ sqrt(sum(pwr(output_resampled))) .* output_resampled;

end


function p = pwr(x)
p = mean(abs(x).^2, 1);
end


function validate_inputs(input, fs, channel)
% Channel version checking
assert(channel.version >= 1.0, ...
    'The minimum version of the channel matrix is 1.0, and you have %.1f.', ...
    channel.version);

% Get constants
T = size(input, 1);
T = T / fs;
T_max = size(channel.h_hat, 3);
T_max = T_max / channel.params.fs_time;

% Check the signal lengths
assert(T < T_max, ...
    'Duration of the input signal, %.2fms, should be no larger than %.2fms.', ...
    T*1e3, T_max*1e3);
end

% [EOF]
