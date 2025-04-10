classdef testReplay < matlab.unittest.TestCase

    properties (TestParameter)
        params = { ...
            struct( ...
            'channel_time', 5, ...
            'coeff', 2, ...
            'fc', 10e3, ...
            'fs_delay', 8e3, ...
            'fs_time', 20, ...
            'has_f_resamp', false, ...
            'has_theta_hat', true, ...
            'n_path', 8, ...
            'R', 4e3, ...
            'Tmp', 10e-3, ...
            'velocity', 1.5), ...
            struct( ...
            'channel_time', 5, ...
            'coeff', 1, ...
            'fc', 15e3, ...
            'fs_delay', 16e3, ...
            'fs_time', 20, ...
            'has_f_resamp', false, ...
            'has_theta_hat', true, ...
            'n_path', 10, ...
            'R', 8e3, ...
            'Tmp', 20e-3, ...
            'velocity', -0.5), ...
            struct( ...
            'channel_time', 5, ...
            'coeff', 1.5, ...
            'fc', 15e3, ...
            'fs_delay', 8e3, ...
            'fs_time', 10, ...
            'has_f_resamp', false, ...
            'has_theta_hat', false, ...
            'n_path', 10, ...
            'R', 4e3, ...
            'Tmp', 20e-3, ...
            'velocity', 0), ...
            struct( ...
            'channel_time', 5, ...
            'coeff', 1.5, ...
            'fc', 15e3, ...
            'fs_delay', 16e3, ...
            'fs_time', 10, ...
            'has_f_resamp', true, ...
            'has_theta_hat', false, ...
            'n_path', 9, ...
            'R', 8e3, ...
            'Tmp', 20e-3, ...
            'velocity', -1), ...
            struct( ...
            'channel_time', 5, ...
            'coeff', 1.5, ...
            'fc', 15e3, ...
            'fs_delay', 10e3, ...
            'fs_time', 10, ...
            'has_f_resamp', true, ...
            'has_theta_hat', true, ...
            'n_path', 8, ...
            'R', 6e3, ...
            'Tmp', 20e-3, ...
            'velocity', -2), ...
            };

    end

    methods (TestMethodSetup)
        function setRandomSeed(~)
            rng(1994);
        end
    end

    methods (Test)

        function testReplayFunction(testCase, params)
            % close all;

            %% Random channel
            path_delay = [0, randsamples(pi/3:params.Tmp*1e3, params.n_path-1)].' / 1e3;
            path_delay = path_delay - min(path_delay);
            path_gain = exp(-path_delay*params.coeff./params.Tmp);
            c_p = path_gain .* exp(-1j*2*pi*(params.fc)*path_delay);

            figure
            h1 = subplot(211);
            stem(path_delay*1e3, path_gain)
            ylabel('Path gain')
            %% Populate the channel matrix

            h_hat = zeros(ceil(params.fs_delay*params.Tmp*1.5), 1, round(params.channel_time*params.fs_time));
            h_hat(round((path_delay + 0.2 * params.Tmp)*params.fs_delay), 1, :) = repmat(c_p, 1, size(h_hat, 3));

            f_resamp = 1/(1 + params.velocity / 1545);
            a = 1 - 1 / f_resamp;
            t = 1:round(params.channel_time*params.fs_delay);
            theta_hat = -a * 2 * pi * params.fc * t / params.fs_delay;

            channel = struct;
            channel.h_hat = h_hat;
            channel.params.fs_delay = params.fs_delay;
            channel.params.fs_time = params.fs_time;
            channel.params.fc = params.fc;
            channel.version = 1.0;

            % If there is additional parameter to resample
            if params.has_theta_hat
                channel.theta_hat = theta_hat;
            end
            if params.has_f_resamp
                channel.f_resamp = f_resamp;
            end
            if params.has_f_resamp && params.has_theta_hat
                channel.theta_hat = zeros(size(channel.theta_hat));
            end

            %% Generate signal
            array_index = 1;
            fs = 48e3;
            data_symbol = randi([0, 1].', 4095, 1) * 2 - 1;
            baseband = resample(data_symbol, fs/params.R, 1);
            passband = real(baseband.*exp(1i*2*pi*params.fc*(0:length(baseband) - 1).'/fs));
            input = [zeros(round(fs/10), 1); passband; zeros(round(fs/10), 1);];

            %% Replay and generate noise
            r = replay(input, fs, array_index, channel);
            % pwelch(r, kaiser(1024, 5), 512, 8192, fs)

            %% Plot the correlation
            h2 = subplot(212);
            v = r .* exp(-2j*pi*params.fc*(0:size(r, 1) - 1).'./fs);
            baseband_resampled = baseband .* exp(-2j*a*pi*params.fc*(0:length(baseband) - 1).'/fs);
            [p, q] = rat(f_resamp);
            baseband_resampled = resample(baseband_resampled, p, q);
            [xcor, lags] = xcorr(v, baseband_resampled);

            xcor(lags <= 0) = [];
            lags(lags <= 0) = [];
            [~, sync] = max(abs(xcor));
            lags = lags - sync;
            xcor = abs(xcor) ./ max(abs(xcor));

            window = (lags >= params.Tmp * -0.2 * fs) .* (lags <= params.Tmp * 1.5 * fs);
            xcor(~window) = [];
            lags(~window) = [];

            [~, locs] = findpeaks(xcor, 'NPeaks', params.n_path, ...
                'MinPeakHeight', min(abs(path_gain))*0.8, ...
                'MinPeakDistance', min(diff(sort(path_delay)))*fs*0.9);

            xaxis = lags.' / fs;
            plot(xaxis*1e3, xcor);
            hold on;
            plot(xaxis(locs)*1e3, xcor(locs), 'x')
            xlabel('Delay [ms]'), ylabel('Xcorr')
            linkaxes([h1, h2], 'xy')

            %% Assertion
            est_delays = xaxis(locs);
            est_gains = xcor(locs);
            [~, index_est] = sort(est_delays);
            [~, index_truth] = sort(path_delay);

            criteria = abs(sum(est_delays(index_est).*est_gains(index_est))- ...
                sum(path_delay(index_truth).*path_gain(index_truth))) < 2e-4 * params.n_path;
            if criteria
                sgtitle('Passed');
            else
                sgtitle('Failed');
            end
            testCase.verifyTrue(criteria == true);
        end

    end
end
function samples = randsamples(population, num)
rand_index = randperm(length(population));
samples = population(rand_index(1:num));
end
