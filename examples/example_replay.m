% Demonstration of the replay and noisegen functions.
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

clc;
clear;
close all;

%% Add the toolbox to the path
addpath('../src');

%% Load channel impulse responses and noise statistics. Refer to README.md for instructions.
channel = load('blue_1.mat');
noise = load('blue_1_noise.mat');

%% Parameters
fs = 48e3; % Sampling rate
fc = 13e3; % Center frequency
R = 4e3; % Symbol rate
M = size(channel.h_hat, 2); % Number of channels
n_repeat = 10; % Number of repeats
array_index = [1, 2, 3]; % Channel index

%% Generate single carrier signals
data_symbol = randi([0, 1].', 1023, 1) * 2 - 1;
baseband = resample(repmat(data_symbol, n_repeat, 1), fs/R, 1);
passband = real(baseband.*exp(1i*2*pi*fc*(0:length(baseband) - 1).'/fs));
input = [zeros(round(fs/10), 1); passband; zeros(round(fs/10), 1);];

%% Replay and generate noise
y = replay(input, fs, array_index, channel);
w = noisegen(size(y), fs);
% w = noisegen(size(y), fs, array_index, noise);

%% Add the noise
r = y + 0.05 * w;

%% Plot the correlation
figure, hold on
legends = cell(length(array_index), 1);
for m = 1:length(array_index)
    v = r(:, m) .* exp(-2j*pi*fc*(0:size(r, 1) - 1).'./fs);
    plot(abs(xcorr(v, resample(data_symbol(1:128), fs/R, 1))));
    legends{m} = sprintf('Receiver %d', array_index(m));
end
xlabel('Samples'), ylabel('Xcorr'), legend(legends)

%% Plot the time domain
t = (0:size(r, 1) - 1) ./ fs;
figure, plot(t, r), legend(legends), xlabel('Time [s]'), ylabel('Received signal')

%% Plot the Welch spectrum
figure, pwelch(r, kaiser(1024, 5), 512, 4096, fs)
xlim([fc - R, fc + R]/1e3)
legend(legends)

%% Remove the toolbox from path
rmpath('../src')

% [EOF]
