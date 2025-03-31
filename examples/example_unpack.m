% Demonstration of the unpack function.
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

%% Load channel impulse responses
channel = load('blue_1.mat');
% channel = rmfield(channel, 'theta_hat');
% channel.f_resamp = 1.0002;

%% Parameters
fs_delay = channel.params.fs_delay;
fs_time = fs_delay * 0.005;
array_index = [1, 3, 5];

%% Unpack the channel
% Buffer fraction for impulse response padding. Increase these values if the
% impulse responses slide out of the window.
buffer_left = 0.1;
buffer_right = 0.1;
unpacked = unpack(fs_time, array_index, channel, buffer_left, buffer_right);

%% Visualize
delay_axis = (0:size(unpacked, 1) - 1) ./ fs_delay;
time_axis = (0:size(unpacked, 3) - 1) ./ fs_time;

figure
imagesc(delay_axis*1e3, time_axis, 20*log10(squeeze(abs(unpacked(:, 1, :))).'), [-30, 0])
xlabel('Delay [ms]')
ylabel('Time [s]')

%% Remove the toolbox from the path
rmpath('../src')

% [EOF]
