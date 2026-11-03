classdef testReplay < matlab.unittest.TestCase
    % TESTREPLAY Unit tests for replay.m
    %
    % Tests {static, mobile} x {theta_hat, phi_hat} across various speed
    % configurations, and f_resamp (bulk Doppler, Watermark) with one factor
    % per hydrophone.
    %
    %   theta_hat: phase tracking only. h_hat drifts with motion.
    %   phi_hat:   phase + delay tracking. h_hat is static.
    %   f_resamp:  testVectorFResamp, testVectorFResampSharedCIR and
    %              testFResampLengthChecked.
    %
    % Other m-files required: replay.m
    % Subfunctions: place_taps, compensate_doppler, randsamples
    % Toolbox required: Signal Processing Toolbox (for resample, findpeaks).
    % MAT-files required: None.
    %
    % See also: replay.m, est.m
    %
    % Author: Zhengnan Li, Claude
    % Email : uwa-channels@ofdm.link
    %
    % License: MIT
    %
    % Revision history:
    %   - Feb. 27, 2026: Initial release.
    %   - Jul. 10, 2026: Added testReplayPower: channel-isolation and
    %                    output-power checks for the array_index API.
    %   - Nov.  3, 2026: Added testVectorFResamp and testFResampLengthChecked
    %                    for one resampling factor per hydrophone, and
    %                    testVectorFResampSharedCIR, which draws them.
    %                    Removed the scalar f_resamp (Watermark) cases.

    properties (Constant)
        c = 1500;
        fs = 48e3;
        start = 1000;
    end

    properties (TestParameter)
        params = { ...
            ...
            ... %% === Static ===
            ...
            struct( ... % Case 1
            'label', 'static_theta', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 20, ...
            'M', 3, 'n_path', 8, 'R', 4e3, 'Tmp', 15e-3, ...
            'coeff', 1.5, 'd', 0.5, 'channel_time', 5, ...
            'tracking', 'theta', ...
            'v_const', 0, 'v_amp', 0, 'n_cycles', 0), ...
            struct( ... % Case 2
            'label', 'static_phi', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 20, ...
            'M', 3, 'n_path', 8, 'R', 4e3, 'Tmp', 15e-3, ...
            'coeff', 1.5, 'd', 0.5, 'channel_time', 5, ...
            'tracking', 'phi', ...
            'v_const', 0, 'v_amp', 0, 'n_cycles', 0), ...
            ...
            ... %% === Low speed AUV (1 m/s) + mild sway ===
            ...
            struct( ... % Case 3
            'label', 'low_speed_theta', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 20, ...
            'M', 3, 'n_path', 8, 'R', 4e3, 'Tmp', 15e-3, ...
            'coeff', 1.5, 'd', 0.5, 'channel_time', 5, ...
            'tracking', 'theta', ...
            'v_const', 1, 'v_amp', 0.2, 'n_cycles', 2), ...
            struct( ... % Case 4
            'label', 'low_speed_phi', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 20, ...
            'M', 3, 'n_path', 8, 'R', 4e3, 'Tmp', 15e-3, ...
            'coeff', 1.5, 'd', 0.5, 'channel_time', 5, ...
            'tracking', 'phi', ...
            'v_const', 1, 'v_amp', 0.2, 'n_cycles', 2), ...
            ...
            ... %% === Moderate speed (3 m/s) + platform sway ===
            ...
            struct( ... % Case 5
            'label', 'moderate_speed_theta', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 100, ...
            'M', 2, 'n_path', 8, 'R', 4e3, 'Tmp', 15e-3, ...
            'coeff', 1.5, 'd', 0.5, 'channel_time', 5, ...
            'tracking', 'theta', ...
            'v_const', 3, 'v_amp', 0.5, 'n_cycles', 3), ...
            struct( ... % Case 6
            'label', 'moderate_speed_phi', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 100, ...
            'M', 2, 'n_path', 8, 'R', 4e3, 'Tmp', 15e-3, ...
            'coeff', 1.5, 'd', 0.5, 'channel_time', 5, ...
            'tracking', 'phi', ...
            'v_const', 3, 'v_amp', 0.5, 'n_cycles', 3), ...
            ...
            ... %% === High speed (5 m/s) + strong sway ===
            ...
            struct( ... % Case 7
            'label', 'high_speed_theta', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 100, ...
            'M', 2, 'n_path', 6, 'R', 4e3, 'Tmp', 15e-3, ...
            'coeff', 1.5, 'd', 0.5, 'channel_time', 5, ...
            'tracking', 'theta', ...
            'v_const', 5, 'v_amp', 1.0, 'n_cycles', 4), ...
            struct( ... % Case 8
            'label', 'high_speed_phi', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 100, ...
            'M', 2, 'n_path', 6, 'R', 4e3, 'Tmp', 15e-3, ...
            'coeff', 1.5, 'd', 0.5, 'channel_time', 5, ...
            'tracking', 'phi', ...
            'v_const', 5, 'v_amp', 1.0, 'n_cycles', 4), ...
            ...
            ... %% === Negative drift (closing) + sway ===
            ...
            struct( ... % Case 9
            'label', 'closing_theta', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 100, ...
            'M', 2, 'n_path', 8, 'R', 4e3, 'Tmp', 15e-3, ...
            'coeff', 1.5, 'd', 0.5, 'channel_time', 5, ...
            'tracking', 'theta', ...
            'v_const', -3, 'v_amp', 0.5, 'n_cycles', 3), ...
            struct( ... % Case 10
            'label', 'closing_phi', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 100, ...
            'M', 2, 'n_path', 8, 'R', 4e3, 'Tmp', 15e-3, ...
            'coeff', 1.5, 'd', 0.5, 'channel_time', 5, ...
            'tracking', 'phi', ...
            'v_const', -3, 'v_amp', 0.5, 'n_cycles', 3), ...
            ...
            ... %% === Pure sway, no drift ===
            ...
            struct( ... % Case 11
            'label', 'sway_only_theta', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 100, ...
            'M', 3, 'n_path', 8, 'R', 4e3, 'Tmp', 15e-3, ...
            'coeff', 1.5, 'd', 0.5, 'channel_time', 5, ...
            'tracking', 'theta', ...
            'v_const', 0, 'v_amp', 1.5, 'n_cycles', 5), ...
            struct( ... % Case 12
            'label', 'sway_only_phi', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 100, ...
            'M', 3, 'n_path', 8, 'R', 4e3, 'Tmp', 15e-3, ...
            'coeff', 1.5, 'd', 0.5, 'channel_time', 5, ...
            'tracking', 'phi', ...
            'v_const', 0, 'v_amp', 1.5, 'n_cycles', 5), ...
            ...
            ... %% === Single element, high speed ===
            ...
            struct( ... % Case 13
            'label', 'single_elem_theta', ...
            'fc', 12e3, 'fs_delay', 8e3, 'fs_time', 100, ...
            'M', 1, 'n_path', 8, 'R', 4e3, 'Tmp', 15e-3, ...
            'coeff', 1.5, 'd', 0.5, 'channel_time', 5, ...
            'tracking', 'theta', ...
            'v_const', 4, 'v_amp', 1.2, 'n_cycles', 3), ...
            struct( ... % Case 14
            'label', 'single_elem_phi', ...
            'fc', 12e3, 'fs_delay', 8e3, 'fs_time', 100, ...
            'M', 1, 'n_path', 8, 'R', 4e3, 'Tmp', 15e-3, ...
            'coeff', 1.5, 'd', 0.5, 'channel_time', 5, ...
            'tracking', 'phi', ...
            'v_const', 4, 'v_amp', 1.2, 'n_cycles', 3), ...
            };
    end

    methods (TestMethodSetup)
        function setRandomSeed(~)
            rng(1994);
        end
    end

    methods (Test)

        function testReplayFunction(testCase, params)

            c = testCase.c;
            fs = testCase.fs;
            start = testCase.start;
            p = params;

            %% 1. Multipath geometry
            path_delay_0 = [0; sort(randsamples( ...
                (pi/3:p.Tmp*1e3) / 1e3, p.n_path - 1)).'];

            incremental_delay = (0:p.M-1) * p.d / c;
            path_delay_0 = path_delay_0 + incremental_delay;
            path_delay_0 = path_delay_0 - min(path_delay_0, [], 'all');

            path_gain = exp(-path_delay_0 * p.coeff ./ p.Tmp);
            c_p = path_gain .* exp(-1j * 2 * pi * p.fc * path_delay_0);

            %% 2. Motion model
            T_ch = p.channel_time;
            N_time = round(T_ch * p.fs_time);
            N_delay = round(T_ch * p.fs_delay);

            t_snapshots = (0:N_time-1) / p.fs_time;
            t_delay = (1:N_delay) / p.fs_delay;

            if p.n_cycles > 0
                omega_osc = 2*pi*p.n_cycles / T_ch;
            end

            dtau_snap = (p.v_const / c) * t_snapshots;
            if p.n_cycles > 0
                dtau_snap = dtau_snap ...
                    - p.v_amp / (c * omega_osc) ...
                    * (cos(omega_osc * t_snapshots) - 1);
            end

            phi = -2*pi*p.fc * (p.v_const/c) .* t_delay;
            phi_sin = zeros(size(t_delay));
            if p.n_cycles > 0
                phi_sin = p.fc * p.v_amp * T_ch / (c * p.n_cycles) ...
                    * (cos(omega_osc * t_delay) - 1);
                phi = phi + phi_sin;
            end

            %% 3. Build h_hat
            L = ceil(p.fs_delay * p.Tmp * 1.5);
            has_motion = (p.v_const ~= 0) || (p.v_amp ~= 0);

            if strcmp(p.tracking, 'theta') && has_motion
                h_hat = zeros(L, p.M, N_time);
                for k = 1:N_time
                    h_hat(:, :, k) = place_taps(L, p.M, ...
                        path_delay_0 + dtau_snap(k), ...
                        c_p, p.Tmp, p.fs_delay);
                end
            else
                h_hat_static = place_taps(L, p.M, ...
                    path_delay_0, c_p, p.Tmp, p.fs_delay);
                h_hat = repmat(h_hat_static, 1, 1, N_time);
            end

            %% 4. Assemble channel struct
            channel = struct();
            channel.h_hat = h_hat;
            channel.params.fs_delay = p.fs_delay;
            channel.params.fs_time = p.fs_time;
            channel.params.fc = p.fc;
            channel.version = 1.0;

            switch p.tracking
                case 'theta'
                    channel.theta_hat = repmat(phi, p.M, 1);
                case 'phi'
                    channel.phi_hat = repmat(phi, p.M, 1);
            end

            %% 5. Transmit signal
            data_symbols = randi([0, 1], 4095, 1) * 2 - 1;
            baseband = resample(data_symbols, fs / p.R, 1);
            passband = real(baseband .* exp(1j*2*pi*p.fc ...
                * (0:length(baseband)-1).'/fs));
            input = [zeros(round(fs/10), 1); passband; zeros(round(fs/10), 1)];

            %% 6. Replay
            r = replay(input, fs, 1:p.M, channel, start);

            %% 7. Phase field for compensation
            switch p.tracking
                case 'theta'
                    phi_field = channel.theta_hat;
                case 'phi'
                    phi_field = channel.phi_hat;
                otherwise
                    phi_field = zeros(p.M, length(t_delay));
            end

            %% 8. Cross-correlation and verification
            baseband_ref = baseband;

            % Replay time window
            [p_rs, q_rs] = rat(p.fs_delay / fs);
            T_baseband = length(resample(baseband, p_rs, q_rs));
            t_replay_start = (start - 1) / p.fs_delay;
            t_replay_end = (start + T_baseband + L) / p.fs_delay;

            figure('Name', p.label);

            % --- h_hat waterfall ---
            subplot(221);
            delay_axis = (0:L-1) / p.fs_delay * 1e3;
            imagesc(delay_axis, t_snapshots, ...
                squeeze(abs(channel.h_hat(:, 1, :))).');
            ylabel('Time [s]'); xlabel('Delay [ms]');
            title('|h_{hat}| (element 1)');
            colorbar;
            hold on;
            yline(t_replay_start, 'r--', 'LineWidth', 1.5);
            yline(t_replay_end, 'r--', 'LineWidth', 1.5);

            % --- Speed from phase ---
            subplot(222);
            dphi = diff(phi_field(1, :));
            dt = 1 / p.fs_delay;
            v_inst = -dphi / (dt * 2 * pi * p.fc) * c;
            t_speed = (0:length(v_inst)-1) / p.fs_delay;
            plot(t_speed, v_inst);
            % A constant speed differs from sample to sample only by rounding
            % error (~1e-12 m/s); a span of at least 1 m/s keeps the axis from
            % zooming in on it.
            v_mid = (max(v_inst) + min(v_inst)) / 2;
            v_half = max(0.55 * (max(v_inst) - min(v_inst)), 0.5);
            ylim(v_mid + [-v_half, v_half]);
            xlabel('Time [s]'); ylabel('Speed [m/s]');
            switch p.tracking
                case 'theta'
                    title('Instantaneous speed (from \theta)');
                case 'phi'
                    title('Instantaneous speed (from \phi)');
                otherwise
                    title('Instantaneous speed (no tracking field)');
            end
            grid on;
            hold on;
            xline(t_replay_start, 'r--', 'LineWidth', 1.5);
            xline(t_replay_end, 'r--', 'LineWidth', 1.5);

            % --- Ground truth + Cross-correlation (merged) ---
            subplot(2, 1, 2); hold on;
            colors = lines(p.M);
            max_xcor = 0;
            sync_idx = zeros(1, p.M);
            criteria = false(1, p.M);

            % Plot ground truth stems first
            for m = 1:p.M
                stem(path_delay_0(:, m)*1e3, path_gain(:, m), ...
                    'Color', colors(m, :), 'MarkerSize', 6);
            end

            for m = 1:p.M
                v_m = r(:, m) .* exp(-2j*pi*p.fc * (0:size(r,1)-1).'/fs);

                if has_motion
                    v_m = compensate_doppler(v_m, fs, p.fc, ...
                        p.fs_delay, phi_field(m, :), start);
                end

                [xcor, lags] = xcorr(v_m, baseband_ref);

                xcor(lags <= 0) = [];
                lags(lags <= 0) = [];

                [~, sync_idx(m)] = max(abs(xcor));

                if m == 1
                    sync_ref = sync_idx(1);
                    max_xcor = max(abs(xcor));
                end
                lags_shifted = lags - sync_ref;
                xcor_norm = abs(xcor) / max_xcor;

                win = (lags_shifted >= -0.2*p.Tmp*fs) ...
                    & (lags_shifted <=  1.5*p.Tmp*fs);
                xcor_win = xcor_norm(win);
                lags_win = lags_shifted(win);

                min_gain = min(abs(path_gain(:, m))) * 0.6;
                min_sep = min(diff(sort(path_delay_0(:, m)))) * fs * 0.7;
                [~, locs] = findpeaks(xcor_win, ...
                    'NPeaks', p.n_path, ...
                    'MinPeakHeight', min_gain, ...
                    'MinPeakDistance', min_sep);

                xaxis = lags_win / fs * 1e3;
                plot(xaxis, xcor_win, '--', 'Color', colors(m, :));
                plot(xaxis(locs), xcor_win(locs), 'x', ...
                    'Color', colors(m, :), 'MarkerSize', 10);

                % Weighted mean delay comparison
                n_found = length(locs);
                est_delays = lags_win(locs) / fs;
                est_gains = xcor_win(locs);
                [~, idx_e] = sort(est_delays);
                [~, idx_t] = sort(path_delay_0(:, m));

                n_found = length(locs);
                if n_found == p.n_path
                    est_delays = sort(lags_win(locs) / fs);
                    true_delays = sort(path_delay_0(:, m));
                    [~, idx_t] = sort(path_delay_0(:, m));
                    true_gains = path_gain(idx_t, m);
                    tol = 2e-4 * p.n_path;
                    criteria(m) = sum(abs(est_delays(:) - true_delays) .* true_gains) < tol;
                end
            end

            xlabel('Delay [ms]'); ylabel('Path gain / |Xcorr|');
            xlim([-0.2*p.Tmp*1e3, 1.5*p.Tmp*1e3]);
            title('Ground truth + Cross-correlation');

            if all(criteria)
                result = 'PASSED';
            else
                result = 'FAILED';
            end
            sgtitle(sprintf('%s: %s', p.label, result), 'Interpreter', 'none');

            testCase.verifyTrue(all(criteria), ...
                sprintf('Peak delay mismatch for case: %s', p.label));
        end

        function testReplayPower(testCase)
            % Two checks for the array_index API:
            %   1. Isolation: replaying through a single hydrophone must equal the
            %      corresponding column of a replay through all hydrophones (same
            %      start). This holds only because output normalization was removed,
            %      so each channel is processed independently of the selected set.
            %   2. Power: after normalizing the channel matrix so the average
            %      per-channel energy is 1 (cf. calib.m), a unit-variance passband
            %      input replays to an output whose per-channel variance is O(1).

            c = testCase.c;
            fs = testCase.fs;
            start = testCase.start;

            fc = 12e3; fs_delay = 10e3; fs_time = 20;
            M = 4; n_path = 8; R = 4e3; Tmp = 15e-3; coeff = 1.5; d = 0.5;
            T_ch = 5;

            %% 1. Multipath geometry (static channel)
            path_delay_0 = [0; sort(randsamples( ...
                (pi/3:Tmp*1e3) / 1e3, n_path - 1)).'];
            incremental_delay = (0:M-1) * d / c;
            path_delay_0 = path_delay_0 + incremental_delay;
            path_delay_0 = path_delay_0 - min(path_delay_0, [], 'all');
            path_gain = exp(-path_delay_0 * coeff ./ Tmp);
            c_p = path_gain .* exp(-1j * 2 * pi * fc * path_delay_0);

            L = ceil(fs_delay * Tmp * 1.5);
            N_time = round(T_ch * fs_time);
            h_hat_static = place_taps(L, M, path_delay_0, c_p, Tmp, fs_delay);
            h_hat = repmat(h_hat_static, 1, 1, N_time);

            %% 2. Normalize so the average per-channel energy is 1 (cf. calib.m)
            T = size(h_hat, 3);
            E_tot = sum(abs(h_hat(:)).^2);
            E_tot_avg = E_tot / T;
            alpha = sqrt(M / E_tot_avg);
            h_hat = alpha * h_hat;

            channel = struct();
            channel.h_hat = h_hat;
            channel.params.fs_delay = fs_delay;
            channel.params.fs_time = fs_time;
            channel.params.fc = fc;
            channel.version = 1.0;

            %% 3. Unit-variance passband probe
            data_symbols = randi([0, 1], 8192, 1) * 2 - 1;
            baseband = resample(data_symbols, fs / R, 1);
            passband = real(baseband .* exp(1j*2*pi*fc*(0:length(baseband)-1).'/fs));
            passband = passband ./ sqrt(var(passband));

            %% 4. Replay one channel, and all channels, with the same start
            target_ch = 3;
            r_single = replay(passband, fs, target_ch, channel, start);
            r_all = replay(passband, fs, 1:M, channel, start);

            %% 5. Check 1: isolation
            n = min(size(r_single, 1), size(r_all, 1));
            ref = r_all(1:n, target_ch);
            err = max(abs(r_single(1:n, 1) - ref)) / max(abs(ref));
            testCase.verifyLessThan(err, 1e-10, ...
                sprintf(['Single-channel replay does not match column %d of the ' ...
                'full replay (relative err = %.2e).'], target_ch, err));

            %% 6. Check 2: per-channel output power is O(1)
            mean_var = mean(var(r_all));
            % fprintf('testReplayPower: per-channel output power mean = %.3f\n', mean_var);
            testCase.verifyGreaterThan(mean_var, 0.93, ...
                sprintf('Mean output power %.3f is unexpectedly small.', mean_var));
            testCase.verifyLessThan(mean_var, 1.07, ...
                sprintf('Mean output power %.3f is unexpectedly large.', mean_var));
        end

        function testLateStartKeepsTheEnd(testCase)
            % A replay late in a drifting channel must not lose either end
            % of the input, whichever way the delay drifts.
            %
            % phi_hat accumulates over the whole record, so late in a channel
            % whose delay drifts steadily its absolute value is large.  replay
            % reinserts the drift relative to the starting sample;
            % reinserting the absolute drift shifted the output by all of it,
            % which cut off the end of the signal when the delay grows
            % (receding) and its beginning when it shrinks (closing).  The
            % delay here drifts by 2 ms per second, so 15 s in it is 30 ms,
            % against 1 ms of change across the 0.5 s input.  A burst at each
            % end of the input must survive.

            fs = 96e3; fc = 13e3; fs_delay = 8e3;
            x = zeros(0.5 * fs, 1);
            b = 0.005 * fs;
            burst = cos(2 * pi * fc * (0:b-1).' / fs);
            x(1:b) = burst;
            x(end-b+1:end) = burst;

            for rate = [2e-3, -2e-3]   % receding, closing
                channel = drift_channel(fc, fs_delay, 400, 64, rate);
                early = replay(x, fs, 1, channel, 1001);
                late = replay(x, fs, 1, channel, 15 * fs_delay + 1);

                ratio = sum(late.^2) / sum(early.^2);
                testCase.verifyGreaterThan(ratio, 0.9, sprintf(['Delay rate ', ...
                    '%+g: a replay 15 s in carries %.2f of the energy of one ', ...
                    'at the start; part of the signal was cut off.'], rate, ratio));
                testCase.verifyLessThan(ratio, 1.1);
                w = 0.05 * fs;
                head = sum(late(1:w).^2) / sum(late.^2);
                tail = sum(late(end-w+1:end).^2) / sum(late.^2);
                testCase.verifyGreaterThan(head, 0.3, sprintf(['Delay rate ', ...
                    '%+g: the burst at the start of the input is missing ', ...
                    '(%.2f).'], rate, head));
                testCase.verifyGreaterThan(tail, 0.3, sprintf(['Delay rate ', ...
                    '%+g: the burst at the end of the input is missing ', ...
                    '(%.2f).'], rate, tail));
            end
        end

        function testOutputLength(testCase)
            % The output is T + L samples at fs_delay, T the input length at
            % fs_delay and L the number of taps, plus, with phi_hat, room for
            % the largest delay the drift reaches in the replay window,
            % measured from the starting sample.  The Python package returns
            % the same length.
            %
            % With a delay growing at a constant rate r, the drift after k
            % samples is r*k samples, so the room needed is ceil(r*(T+L-1));
            % a shrinking delay needs none, and so do theta_hat and no
            % tracking.

            fs = 96e3; fc = 13e3; fs_delay = 8e3; L = 64;
            x = randn(0.5 * fs, 1);
            T = ceil(numel(x) * fs_delay / fs);
            up = fs / fs_delay;
            rate = 2e-3;

            channel = drift_channel(fc, fs_delay, 400, L, rate);
            testCase.verifyEqual(numel(replay(x, fs, 1, channel, 1001)), ...
                (T + L + ceil(rate * (T + L - 1))) * up, ...
                'phi_hat, growing delay: output length');

            channel = drift_channel(fc, fs_delay, 400, L, -rate);
            testCase.verifyEqual(numel(replay(x, fs, 1, channel, 1001)), ...
                (T + L) * up, 'phi_hat, shrinking delay: output length');

            channel.theta_hat = channel.phi_hat;
            channel = rmfield(channel, 'phi_hat');
            testCase.verifyEqual(numel(replay(x, fs, 1, channel, 1001)), ...
                (T + L) * up, 'theta_hat: output length');

            channel = rmfield(channel, 'theta_hat');
            testCase.verifyEqual(numel(replay(x, fs, 1, channel, 1001)), ...
                (T + L) * up, 'no tracking: output length');
        end

        function testVectorFResamp(testCase)
            % f_resamp may hold one factor per hydrophone.  Each column of the
            % output must then equal a replay through that hydrophone alone
            % with its factor as a scalar, zero-padded to the longest column.
            % A vector of equal factors must give the same output as the
            % scalar, and array_index must pick the factors of the
            % hydrophones it selects, in its order.

            fs = 48e3; start = testCase.start;
            fc = 12e3; fs_delay = 10e3; fs_time = 20; L = 64; N_time = 100;
            M = 3;
            h_hat = zeros(L, M, N_time);
            h_hat(10, 1, :) = 1;
            h_hat([12, 30], 2, :) = repmat([0.8; 0.4j], 1, 1, N_time);
            h_hat([8, 21, 40], 3, :) = repmat([1; -0.5; 0.3j], 1, 1, N_time);
            channel = struct();
            channel.h_hat = h_hat;
            channel.version = 1.0;
            channel.params = struct('fs_delay', fs_delay, 'fs_time', fs_time, 'fc', fc);
            t = (0:N_time/fs_time*fs_delay - 1) / fs_delay;
            channel.phi_hat = repmat(-2 * pi * fc * 1e-4 * t, M, 1);

            x = cos(2 * pi * fc * (0:0.5*fs-1).' / fs) .* hann(0.5*fs);
            f_resamp = 1 ./ (1 + [2e-3; -1e-3; 4e-3]);

            % Each column against a single-hydrophone, scalar-factor replay.
            channel.f_resamp = f_resamp;
            for index = {1:M, [3, 1]}
                array_index = index{1};
                r = replay(x, fs, array_index, channel, start);
                lengths = zeros(1, numel(array_index));
                for k = 1:numel(array_index)
                    m = array_index(k);
                    single = channel;
                    single.h_hat = channel.h_hat(:, m, :);
                    single.phi_hat = channel.phi_hat(m, :);
                    single.f_resamp = f_resamp(m);
                    r_m = replay(x, fs, 1, single, start);
                    lengths(k) = numel(r_m);
                    n = lengths(k);
                    err = max(abs(r(1:n, k) - r_m)) / max(abs(r_m));
                    testCase.verifyLessThan(err, 1e-10, sprintf(['Column %d ', ...
                        '(hydrophone %d) of replay(..., %s, ...) differs from ', ...
                        'hydrophone %d replayed alone with f_resamp = %.6f ', ...
                        '(relative err = %.2e).'], k, m, mat2str(array_index), ...
                        m, f_resamp(m), err));
                    testCase.verifyEqual(r(n+1:end, k), zeros(size(r, 1) - n, 1), ...
                        sprintf('Column %d is not zero-padded past its end.', k));
                end
                testCase.verifyEqual(size(r, 1), max(lengths), ...
                    'The output is not as long as its longest column.');
                testCase.verifyGreaterThan(numel(unique(lengths)), 1, ...
                    'The factors are too close to change the column lengths.');
            end

            % Equal factors as a row vector give the scalar result.
            channel.f_resamp = f_resamp(1);
            r_scalar = replay(x, fs, 1:M, channel, start);
            channel.f_resamp = repmat(f_resamp(1), 1, M);
            r_vector = replay(x, fs, 1:M, channel, start);
            testCase.verifyEqual(size(r_vector), size(r_scalar));
            testCase.verifyLessThan(max(abs(r_vector - r_scalar), [], 'all') ...
                / max(abs(r_scalar), [], 'all'), 1e-12, ...
                'A vector of equal factors differs from the scalar factor.');
        end

        function testVectorFResampSharedCIR(testCase)
            % Three hydrophones share one static CIR and differ only in
            % f_resamp, so any difference between their outputs comes from
            % f_resamp alone.  replay resamples column m by f_resamp(m), which
            % stretches it in time by that factor: a column is f_resamp(m)
            % times as long as without it, and two events D seconds apart
            % come out f_resamp(m)*D apart.
            %
            % The input is two short bursts at fc, D seconds apart.  The test
            % checks each column's length and burst spacing against the
            % factor, and the figure shows them: on the left the whole column,
            % with its end marked and its length labelled; on the right the
            % second burst, with where it falls without f_resamp (dashed) and
            % how far f_resamp moves it.

            fs = 48e3; start = testCase.start;
            fc = 12e3; fs_delay = 8e3; fs_time = 20; L = 64; N_time = 100;
            M = 3;
            taps = zeros(L, 1);
            taps([10, 18, 30]) = [1; 0.5 - 0.3j; 0.25j];
            channel = struct();
            channel.h_hat = repmat(taps, [1, M, N_time]);
            channel.version = 1.0;
            channel.params = struct('fs_delay', fs_delay, 'fs_time', fs_time, 'fc', fc);

            % Two 2-ms bursts, D = 0.4 s apart, in 0.5 s of input.
            b = 2e-3 * fs;
            burst = cos(2 * pi * fc * (0:b-1).' / fs) .* hann(b);
            x = zeros(0.5 * fs, 1);
            n1 = 0.05 * fs; D = 0.4;
            x(n1 + (1:b)) = burst;
            x(n1 + D * fs + (1:b)) = burst;

            % Without f_resamp, as the reference.
            r0 = replay(x, fs, 1, channel, start);
            q0 = burst_peaks(r0);
            len0 = numel(r0);

            excess = [1e-2; -1e-2; -3e-2];               % 1/f_resamp - 1
            f_resamp = 1 ./ (1 + excess);
            channel.f_resamp = f_resamp;
            r = replay(x, fs, 1:M, channel, start);

            figure('Name', 'vector_f_resamp_replay');
            lengths = zeros(M, 1);
            for m = 1:M
                single = channel;
                single.h_hat = channel.h_hat(:, m, :);
                single.f_resamp = f_resamp(m);
                lengths(m) = numel(replay(x, fs, 1, single, start));
                q = burst_peaks(r(1:lengths(m), m));

                testCase.verifyLessThan(abs(lengths(m) - len0 * f_resamp(m)), 2, ...
                    sprintf(['Hydrophone %d: the column is %d samples long, ', ...
                    'not f_resamp * %d = %.1f.'], m, lengths(m), len0, ...
                    len0 * f_resamp(m)));
                testCase.verifyLessThan(abs(diff(q) - diff(q0) * f_resamp(m)), 2, ...
                    sprintf(['Hydrophone %d: the bursts are %d samples apart, ', ...
                    'not f_resamp * %d = %.1f.'], m, diff(q), diff(q0), ...
                    diff(q0) * f_resamp(m)));

                env = abs(hilbert(r(:, m)));
                env = env / max(env);
                t = (0:numel(env) - 1).' / fs;

                % Whole column, its end marked and its length labelled.
                subplot(M, 2, 2*m - 1);
                plot(t, env); hold on;
                xline(lengths(m) / fs, 'r-', 'LineWidth', 1.5);
                xline(len0 / fs, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
                xlim([0, size(r, 1) / fs * 1.02]); ylim([0, 1.1]);
                ylabel('|r|'); grid on;
                title(sprintf(['Hydrophone %d: f_{resamp} = 1/(1%+.0e), ', ...
                    'length %d (%d without)'], m, excess(m), lengths(m), len0));
                if m == M, xlabel('Time [s]'); end

                % Second burst, without f_resamp (dashed) and with it.
                subplot(M, 2, 2*m);
                plot(t * 1e3, env); hold on;
                xline(q0(2) / fs * 1e3, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
                xline(q(2) / fs * 1e3, 'r-');
                shift = (q(2) - q0(2)) / fs * 1e3;
                plot([q0(2), q(2)] / fs * 1e3, [1.05, 1.05], 'r-', 'LineWidth', 2);
                text((q0(2) + q(2)) / 2 / fs * 1e3, 0.85, ...
                    sprintf('%+.2f ms', shift), 'Color', 'r', ...
                    'HorizontalAlignment', 'center', 'FontWeight', 'bold');
                xlim(q0(2) / fs * 1e3 + [-8, 20]); ylim([0, 1.15]);
                grid on;
                title(sprintf('Second burst moves by %+.2f ms', shift));
                if m == M, xlabel('Time [ms]'); end
            end
            sgtitle(['Same CIR, different f_{resamp}: each column is stretched ', ...
                'by its own factor (dashed: without f_{resamp})']);

            testCase.verifyEqual(numel(unique(lengths)), M, ...
                'The columns should all have different lengths.');
            testCase.verifyEqual(size(r, 1), max(lengths), ...
                'The output is not as long as its longest column.');
        end

        function testFResampLengthChecked(testCase)
            % A vector f_resamp must have one factor per hydrophone in the
            % file, whatever array_index selects.
            channel = drift_channel(13e3, 8e3, 40, 64, 1e-4);
            channel.h_hat = repmat(channel.h_hat, 1, 3, 1);
            channel.phi_hat = repmat(channel.phi_hat, 3, 1);
            channel.f_resamp = [1, 1];
            x = randn(0.1 * 96e3, 1);
            testCase.verifyError(@() replay(x, 96e3, 1, channel, 1001), ...
                'replay:invalidFResamp');
            channel.f_resamp = [1, -1, 1];
            testCase.verifyError(@() replay(x, 96e3, 1, channel, 1001), ...
                'replay:invalidFResamp');
        end

        function testRandomStartStaysInRange(testCase)
            % Without a start, replay draws one at random.  With the input
            % nearly as long as the channel allows, the valid starts are a
            % handful of samples, so a draw past the last one reads beyond
            % phi_hat or h_hat and errors, and one that fell short of the
            % first would too.

            fs = 96e3; fc = 13e3; fs_delay = 8e3; fs_time = 20; L = 64;
            N_time = 40;                                % 2 s of channel
            channel = drift_channel(fc, fs_delay, N_time, L, 2e-3);
            T_max = N_time / fs_time * fs_delay;
            T = T_max - L - 5;                          % 5 valid starts
            x = randn(T * fs / fs_delay, 1);
            for k = 1:200
                y = replay(x, fs, 1, channel);
                testCase.verifyTrue(all(isfinite(y)), 'non-finite output');
            end
        end

    end
end


%% Helpers

function peaks = burst_peaks(r)
% Sample indices of the two burst peaks in R, one in each half of it.
env = abs(hilbert(r));
half = round(numel(env) / 2);
[~, p1] = max(env(1:half));
[~, p2] = max(env(half+1:end));
peaks = [p1, p2 + half];
end

function channel = drift_channel(fc, fs_delay, N_time, L, rate)
% A single-tap channel whose delay drifts at RATE seconds per second, held
% as phi_hat = -2*pi*fc*rate*t, with fs_time = 20 Hz.
fs_time = 20;
h_hat = zeros(L, 1, N_time);
h_hat(L/2, 1, :) = 1;
channel = struct();
channel.h_hat = h_hat;
channel.version = 1.0;
channel.params = struct('fs_delay', fs_delay, 'fs_time', fs_time, 'fc', fc);
t = (0:N_time/fs_time*fs_delay - 1) / fs_delay;
channel.phi_hat = -2 * pi * fc * rate * t;
end

function h = place_taps(L, M, delays, gains, Tmp, fs_delay)
h = zeros(L, M);
subs = round((delays + 0.2*Tmp) * fs_delay);
valid = (subs >= 1) & (subs <= L);
col_offset = (0:M-1) * L;
for m = 1:M
    v = valid(:, m);
    h(subs(v, m) + col_offset(m)) = gains(v, m);
end
end


function v_out = compensate_doppler(v_in, fs, fc, fs_delay, phi_field, start)
N = length(v_in);
t_rx = (0:N-1).' / fs;

t_phi = (0:length(phi_field)-1) / fs_delay;
t_abs = t_rx + (start - 1) / fs_delay;

phi_rx = interp1(t_phi, phi_field, t_abs, 'spline');
phi_rx = phi_rx - phi_rx(1);

dtau = -phi_rx / (2*pi*fc);

v_comp = v_in .* exp(-1j * phi_rx);
v_out = interp1(t_rx, v_comp, t_rx + dtau, 'spline', 0);
end


function samples = randsamples(population, num)
rand_index = randperm(length(population));
samples = population(rand_index(1:num));
end
