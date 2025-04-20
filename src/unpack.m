function unpacked_channel = unpack(fs, array_index, channel, varargin)
% UNPACK Unpack the channel coefficients.
% Uncompress the channel coefficients.
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
%    see example_unpack.m.
%
% Other m-files required: None
% Subfunctions: None
% Toolbox required: Signal Processing Toolbox (R) (resample function).
% MAT-files required: channel matfile.
%
% See also:
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

%% Unpacking variables
params = channel.params;
h_hat = channel.h_hat(:, array_index, :);

%% Parameters
fs_delay = params.fs_delay; % Sample rate in delay
fs_time = params.fs_time; % Sample rate in time
fs_time_desired = fs;
fc = params.fc; % Center frequency
K = size(h_hat, 1); % Delay axis
M = size(h_hat, 2); % Number of array elements
T = size(h_hat, 3); % Time axis

if isfield(channel, 'theta_hat')
    theta_hat = channel.theta_hat(array_index, :);
end
if isfield(channel, 'f_resamp') % If we have additional parameter to resample
    theta_hat = repmat((1 / channel.f_resamp - 1) * 2 * pi * fc * (1:ceil(T*fs_delay/fs_time)) / fs_delay, M, 1);
end

%% Allocate some buffer
if nargin == 3
    buffer_left = 0.1; % Ratio of K, allocated to the left
    buffer_right = 0.1; % Ratio of K, allocated to the right
else
    buffer_left = varargin{1};
    buffer_right = varargin{2};
end
h_hat = [zeros(ceil(K*buffer_left), M, T); h_hat; zeros(ceil(K*buffer_right), M, T)];
K = size(h_hat, 1);

%% Sample rate conversion
[p1, q1] = rat(fs_time_desired/fs_time);
[p2, q2] = rat(fs_time_desired/fs_delay);
delays = (0:K - 1) ./ fs_delay;

%% Unpack the channel for every array element
unpacked_channel = zeros(size(h_hat, 1), length(array_index), ceil(size(h_hat, 3)*p1/q1));
for m = 1:length(array_index)
    h_hat_m = squeeze(h_hat(:, m, :));
    if isfield(channel, 'theta_hat') || isfield(channel, 'f_resamp')
        theta_hat_resampled = resample(theta_hat(m, :), p2, q2);
        drift = theta_hat_resampled ./ (2 * pi * fc);
        unpacked_channel(:, m, :) = resample(h_hat_m, p1, q1, 'Dimension', 2) .* exp(1j*theta_hat_resampled);
        for t = 1:size(unpacked_channel, 3)
            unpacked_channel(:, m, t) = interp1(delays, squeeze(unpacked_channel(:, m, t)), delays+drift(t), 'spline');
        end
    else
        unpacked_channel(:, m, :) = resample(h_hat_m, p1, q1, 'Dimension', 2);
    end
end
unpacked_channel = unpacked_channel ./ max(abs(unpacked_channel), [], 'all');
end

% [EOF]
