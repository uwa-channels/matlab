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

%% Load channel impulse responses
channel = load('blue_1.mat');
% channel = rmfield(channel, 'theta_hat');
% channel.f_resamp = 1.0002;

%% Parameters
fs_delay = channel.params.fs_delay;
fs_time = fs_delay * 0.005;
array_index = [1, 3, 5];

%% Unpack the channel
% The delay axis is padded to hold the delay drift.  To pad by a fraction of
% the taps instead, pass it: unpack(fs_time, array_index, channel, 0.1, 0.1).
unpacked = unpack(fs_time, array_index, channel);

%% Visualize
delay_axis = (0:size(unpacked, 1) - 1) ./ fs_delay;
time_axis = (0:size(unpacked, 3) - 1) ./ fs_time;

figure
imagesc(delay_axis*1e3, time_axis, 20*log10(squeeze(abs(unpacked(:, 1, :))).'), [-30, 0])
xlabel('Delay [ms]')
ylabel('Time [s]')

% [EOF]
