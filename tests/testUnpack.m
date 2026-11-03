classdef testUnpack < matlab.unittest.TestCase
    % TESTUNPACK Unit tests for unpack.m
    %
    % Tests {static, mobile} x {none, theta_hat, phi_hat}, and f_resamp
    % (bulk Doppler, Watermark) with one factor per hydrophone.
    %
    %   theta_hat: phase tracking only. h_hat drifts with motion.
    %     Unpack re-inserts phase only (no delay interpolation).
    %   phi_hat: phase + delay tracking. h_hat is static.
    %     Unpack re-inserts both phase and delay drift.
    %   f_resamp: bulk Doppler. Adds linear phase ramp + delay drift.
    %     testFResampRampOrigin, testFResampRampRate, testVectorFResamp.
    %
    % Other m-files required: unpack.m
    % Subfunctions: place_taps, randsamples
    % Toolbox required: Signal Processing Toolbox (for resample).
    % MAT-files required: None.
    %
    % See also: unpack.m, replay.m
    %
    % Author: Zhengnan Li, Claude
    % Email : uwa-channels@ofdm.link
    %
    % License: MIT
    %
    % Revision history:
    %   - Feb. 27, 2026: Initial release.
    %   - Sep. 16, 2026: Added testFResampRampOrigin and testFResampRampRate.
    %                    The parameterized test above draws pictures and
    %                    asserts nothing, so the f_resamp cases in it passed
    %                    while the ramp started one sample late.
    %   - Nov.  3, 2026: Added testVectorFResamp for one resampling factor
    %                    per hydrophone, and removed the scalar f_resamp
    %                    (Watermark) cases from testUnpackFunction.
    %                    Added testInvalidInputs for validate_inputs.

    properties (Constant)
        c = 1500;
    end

    properties (TestParameter)
        params = { ...
            ...
            ... %% === Static ===
            ...
            struct( ... % Case 1: Static, no tracking
            'label', 'static_none', ...
            'fc', 10e3, 'fs_delay', 8e3, 'fs_time', 20, ...
            'n_path', 8, 'Tmp', 10e-3, 'coeff', 1, ...
            'channel_time', 5, ...
            'tracking', 'none', ...
            'v_const', 0, 'v_amp', 0, 'n_cycles', 0), ...
            ...
            ... %% === Constant drift, theta_hat ===
            ...
            struct( ... % Case 2: theta_hat, moderate speed
            'label', 'drift_theta', ...
            'fc', 10e3, 'fs_delay', 8e3, 'fs_time', 100, ...
            'n_path', 8, 'Tmp', 10e-3, 'coeff', 1, ...
            'channel_time', 5, ...
            'tracking', 'theta', ...
            'v_const', -2, 'v_amp', 0, 'n_cycles', 0), ...
            ...
            struct( ... % Case 3: theta_hat, high speed
            'label', 'drift_theta_fast', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 100, ...
            'n_path', 8, 'Tmp', 15e-3, 'coeff', 1.5, ...
            'channel_time', 5, ...
            'tracking', 'theta', ...
            'v_const', 4, 'v_amp', 0, 'n_cycles', 0), ...
            ...
            ... %% === Constant drift, phi_hat ===
            ...
            struct( ... % Case 4: phi_hat, moderate speed
            'label', 'drift_phi', ...
            'fc', 10e3, 'fs_delay', 8e3, 'fs_time', 20, ...
            'n_path', 8, 'Tmp', 10e-3, 'coeff', 1, ...
            'channel_time', 5, ...
            'tracking', 'phi', ...
            'v_const', -2, 'v_amp', 0, 'n_cycles', 0), ...
            ...
            struct( ... % Case 5: phi_hat, high speed
            'label', 'drift_phi_fast', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 20, ...
            'n_path', 8, 'Tmp', 15e-3, 'coeff', 1.5, ...
            'channel_time', 5, ...
            'tracking', 'phi', ...
            'v_const', 4, 'v_amp', 0, 'n_cycles', 0), ...
            ...
            ... %% === Sinusoidal motion ===
            ...
            struct( ... % Case 6: sway only, theta_hat
            'label', 'sway_theta', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 100, ...
            'n_path', 8, 'Tmp', 15e-3, 'coeff', 1.5, ...
            'channel_time', 5, ...
            'tracking', 'theta', ...
            'v_const', 0, 'v_amp', 1.5, 'n_cycles', 4), ...
            ...
            struct( ... % Case 7: sway only, phi_hat
            'label', 'sway_phi', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 20, ...
            'n_path', 8, 'Tmp', 15e-3, 'coeff', 1.5, ...
            'channel_time', 5, ...
            'tracking', 'phi', ...
            'v_const', 0, 'v_amp', 1.5, 'n_cycles', 4), ...
            ...
            ... %% === Combined: AUV + sway ===
            ...
            struct( ... % Case 8: combined, theta_hat
            'label', 'combined_theta', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 100, ...
            'n_path', 8, 'Tmp', 15e-3, 'coeff', 1.5, ...
            'channel_time', 5, ...
            'tracking', 'theta', ...
            'v_const', 3, 'v_amp', 0.5, 'n_cycles', 3), ...
            ...
            struct( ... % Case 9: combined, phi_hat
            'label', 'combined_phi', ...
            'fc', 12e3, 'fs_delay', 10e3, 'fs_time', 20, ...
            'n_path', 8, 'Tmp', 15e-3, 'coeff', 1.5, ...
            'channel_time', 5, ...
            'tracking', 'phi', ...
            'v_const', 3, 'v_amp', 0.5, 'n_cycles', 3), ...
            };
    end

    methods (TestMethodSetup)
        function setRandomSeed(~)
            rng(1994);
        end
    end

    methods (Test)

        function testUnpackFunction(testCase, params)

            c = testCase.c;
            p = params;

            %% 1. Multipath geometry (single element)
            path_delay_0 = [0; sort(randsamples( ...
                (pi/3:p.Tmp*1e3) / 1e3, p.n_path - 1)).'];
            path_delay_0 = path_delay_0 - min(path_delay_0);

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

            % Cumulative delay at snapshot times
            dtau_snap = (p.v_const / c) * t_snapshots;
            if p.n_cycles > 0
                dtau_snap = dtau_snap ...
                    - p.v_amp / (c * omega_osc) ...
                    * (cos(omega_osc * t_snapshots) - 1);
            end

            % Cumulative phase at fs_delay rate
            %   phi(t) = -2*pi*fc * integral_0^t v(s)/c ds
            phi_const = -2*pi*p.fc * (p.v_const/c) .* t_delay;
            if p.n_cycles > 0
                phi_sin = p.fc * p.v_amp * T_ch / (c * p.n_cycles) ...
                    * (cos(omega_osc * t_delay) - 1);
            else
                phi_sin = zeros(size(t_delay));
            end
            phi_full = phi_const + phi_sin;

            %% 3. Build h_hat
            L = ceil(p.fs_delay * p.Tmp * 1.5);
            has_motion = (p.v_const ~= 0) || (p.v_amp ~= 0);

            if strcmp(p.tracking, 'theta') && has_motion
                % theta_hat: delays drift in h_hat
                h_hat = zeros(L, 1, N_time);
                for k = 1:N_time
                    h_hat(:, 1, k) = place_taps(L, 1, ...
                        path_delay_0 + dtau_snap(k), ...
                        c_p, p.Tmp, p.fs_delay);
                end
            else
                % phi_hat / none / static: h_hat is fixed
                h_hat_static = place_taps(L, 1, ...
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
                    channel.theta_hat = phi_full;
                case 'phi'
                    channel.phi_hat = phi_full;
                case 'none'
                    % No tracking fields
            end

            %% 5. Unpack
            fs_time_out = p.fs_delay * 0.01;
            array_index = 1;
            unpacked = unpack(fs_time_out, array_index, channel, 0.3, 0.3);

            %% 6. Plot
            delay_axis = (0:size(unpacked, 1)-1) ./ p.fs_delay;
            time_axis = (0:size(unpacked, 3)-1) ./ fs_time_out;

            figure('Name', p.label);

            % --- Unpacked channel ---
            subplot(131);
            imagesc(delay_axis*1e3, time_axis, ...
                20*log10(squeeze(abs(unpacked(:, 1, :))).'), [-30, 0]);
            xlabel('Delay [ms]'); ylabel('Time [s]');
            title('Unpacked');
            colorbar;

            % --- h_hat (input) ---
            subplot(132);
            imagesc(delay_axis(1:L)*1e3, t_snapshots, ...
                20*log10(squeeze(abs(h_hat(:, 1, :))).'), [-30, 0]);
            xlabel('Delay [ms]'); ylabel('Time [s]');
            title('h_{hat} (input)');
            colorbar;

            % --- Speed from phase ---
            subplot(133);
            if strcmp(p.tracking, 'theta')
                phi_plot = channel.theta_hat;
            elseif strcmp(p.tracking, 'phi')
                phi_plot = channel.phi_hat;
            else
                phi_plot = zeros(1, N_delay);
            end
            dphi = diff(phi_plot(1, :));
            v_inst = -dphi / (1/p.fs_delay * 2 * pi * p.fc) * c;
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
                    title('Speed (from \theta)');
                case 'phi'
                    title('Speed (from \phi)');
                otherwise
                    title('Speed (no tracking field)');
            end
            grid on;

            sgtitle(p.label, 'Interpreter', 'none');
        end

        function testFResampRampOrigin(testCase)
            % The f_resamp phase ramp must start at zero.
            %
            % unpack builds the ramp on the grid t_orig = (0:N_phi-1)/fs_delay, so
            % sample n of the trajectory is at (n-1)/fs_delay and the ramp is zero
            % at t = 0.  Building it on (1:N_phi) instead puts a constant phase
            % rotation, and through phase_drift a constant delay offset, on every
            % unpacked tap -- which is what this checks for.
            %
            % The assertion compares the same channel unpacked with and without
            % f_resamp at the first output sample.  Both runs share one h_resampled,
            % so resampling transients cancel exactly and the only thing left is the
            % phase and drift that f_resamp contributes, which at t = 0 is nothing.

            fc = 15e3; fs_delay = 16e3; fs_time = 100; fs_out = 160;
            L = 128; N_time = 200;
            excess = 1e-3;                       % 1/f_resamp - 1
            f_resamp = 1 / (1 + excess);

            % A static channel: a handful of taps, identical in every snapshot.
            taps = zeros(L, 1);
            taps([5, 19, 44, 77]) = [1; 0.6 - 0.3j; -0.4 + 0.2j; 0.15j];
            h_hat = repmat(taps, [1, 1, N_time]);

            channel = struct();
            channel.h_hat = h_hat;
            channel.version = 1.0;
            channel.params = struct('fs_delay', fs_delay, 'fs_time', fs_time, 'fc', fc);

            u0 = unpack(fs_out, 1, channel, 0.3, 0.3);
            channel.f_resamp = f_resamp;
            u1 = unpack(fs_out, 1, channel, 0.3, 0.3);

            first0 = squeeze(u0(:, 1, 1));
            first1 = squeeze(u1(:, 1, 1));
            scale = norm(first0);
            testCase.verifyGreaterThan(scale, 0, ...
                'The test channel unpacked to nothing; the fixture is wrong.');

            % What the wrong origin would have done, as a yardstick: a rotation by
            % one sample of the ramp.  The tolerance sits three orders of magnitude
            % below it, and well above floating point.
            wrong = abs(exp(1j * excess * 2 * pi * fc / fs_delay) - 1);
            testCase.verifyLessThan(norm(first1 - first0) / scale, wrong / 1e3, ...
                sprintf(['The f_resamp ramp does not start at zero: the first ', ...
                'output sample moves by %.2e of its norm when f_resamp is added, ', ...
                'against %.2e for a ramp built one sample in.'], ...
                norm(first1 - first0) / scale, wrong));
        end

        function testFResampRampRate(testCase)
            % ... and it must then advance at (1/f_resamp - 1) * 2*pi*fc per second.
            %
            % Measured on the sum over the delay axis, which is the zero-frequency
            % point of the delay response and so is unchanged by the delay drift
            % f_resamp also applies.  excess is kept small enough that the drift
            % stays well inside one sample over the record, so nothing leaves the
            % delay window.

            fc = 15e3; fs_delay = 16e3; fs_time = 100; fs_out = 160;
            L = 128; N_time = 500;
            excess = 2e-6;
            f_resamp = 1 / (1 + excess);

            taps = zeros(L, 1);
            taps([21, 40, 63]) = [1; 0.5 - 0.2j; 0.25 + 0.1j];
            h_hat = repmat(taps, [1, 1, N_time]);

            channel = struct();
            channel.h_hat = h_hat;
            channel.version = 1.0;
            channel.params = struct('fs_delay', fs_delay, 'fs_time', fs_time, 'fc', fc);

            u0 = unpack(fs_out, 1, channel, 0.3, 0.3);
            channel.f_resamp = f_resamp;
            u1 = unpack(fs_out, 1, channel, 0.3, 0.3);

            g0 = squeeze(sum(u0(:, 1, :), 1));
            g1 = squeeze(sum(u1(:, 1, :), 1));
            keep = abs(g0) > 0.1 * max(abs(g0));       % skip the resample ramp-up
            t = (0:numel(g0) - 1).' / fs_out;

            measured = angle(g1(keep) ./ g0(keep));
            predicted = excess * 2 * pi * fc * t(keep);
            residual = angle(exp(1j * (measured - predicted)));   % wrap-safe

            testCase.verifyLessThan(max(abs(residual)), 1e-3, ...
                sprintf(['The f_resamp phase ramp advances at the wrong rate: ', ...
                'worst residual %.2e rad against a total ramp of %.3f rad.'], ...
                max(abs(residual)), predicted(end)));
        end

        function testVectorFResamp(testCase)
            % f_resamp may hold one factor per hydrophone.  Here all three
            % hydrophones share one static CIR and differ only in f_resamp,
            % so any difference between them in the unpacked channel comes
            % from f_resamp alone.  Each hydrophone must unpack as it does
            % alone with its factor as a scalar, array_index must pick the
            % factors of the hydrophones it selects, and a vector of the
            % wrong length must be rejected.
            %
            % The figure shows the three unpacked channels.  A factor
            % f_resamp = 1/(1 + e) moves every arrival by -e*t, so over the
            % record of T seconds the CIR slides by -e*T: earlier for e > 0,
            % later for e < 0.  Each panel marks that slide, its length in ms,
            % and the track the first arrival is predicted to follow.

            fc = 15e3; fs_delay = 16e3; fs_time = 100; fs_out = 160;
            L = 128; N_time = 200; M = 3;
            buf = 0.3;
            taps = zeros(L, 1);
            taps([40, 55, 70]) = [1; 0.6 - 0.3j; 0.3 + 0.2j];
            h_hat = repmat(taps, [1, M, N_time]);

            channel = struct();
            channel.h_hat = h_hat;
            channel.version = 1.0;
            channel.params = struct('fs_delay', fs_delay, 'fs_time', fs_time, 'fc', fc);
            excess = [5e-4; -5e-4; -1.5e-3];             % 1/f_resamp - 1
            f_resamp = 1 ./ (1 + excess);
            channel.f_resamp = f_resamp;

            for index = {1:M, [3, 1]}
                array_index = index{1};
                u = unpack(fs_out, array_index, channel, buf, buf);
                for k = 1:numel(array_index)
                    m = array_index(k);
                    single = channel;
                    single.h_hat = channel.h_hat(:, m, :);
                    single.f_resamp = f_resamp(m);
                    u_m = unpack(fs_out, 1, single, buf, buf);
                    err = max(abs(u(:, k, :) - u_m), [], 'all') ...
                        / max(abs(u_m), [], 'all');
                    testCase.verifyLessThan(err, 1e-12, sprintf(['Hydrophone %d ', ...
                        'in unpack(..., %s, ...) differs from it unpacked alone ', ...
                        'with f_resamp = %.6f (relative err = %.2e).'], ...
                        m, mat2str(array_index), f_resamp(m), err));
                end
            end

            % Same CIR, so with a shared factor the hydrophones would be
            % identical; with their own factors they must all differ.
            u = unpack(fs_out, 1:M, channel, buf, buf);
            for pair = [1 2; 1 3; 2 3].'
                d = max(abs(u(:, pair(1), :) - u(:, pair(2), :)), [], 'all');
                testCase.verifyGreaterThan(d, 0.1, sprintf(['Hydrophones %d ', ...
                    'and %d unpack alike although their f_resamp differ.'], pair));
            end

            % Each hydrophone's first arrival must end the record where its
            % factor puts it, to within a sample.
            n_left = ceil(L * buf);
            delay_axis = ((0:size(u, 1) - 1) - n_left) / fs_delay;
            time_axis = (0:size(u, 3) - 1) / fs_out;
            T = time_axis(end);
            tau0 = (40 - 1) / fs_delay;
            for m = 1:M
                last = abs(u(:, m, end));
                window = abs(delay_axis - (tau0 - excess(m) * T)) < 4 / fs_delay;
                [~, i] = max(last .* window(:));
                testCase.verifyLessThan(abs(delay_axis(i) - (tau0 - excess(m) * T)), ...
                    1 / fs_delay, sprintf(['Hydrophone %d: the first arrival ends ', ...
                    'at %.3f ms, not at %.3f ms.'], m, delay_axis(i) * 1e3, ...
                    (tau0 - excess(m) * T) * 1e3));
            end

            figure('Name', 'vector_f_resamp');
            for m = 1:M
                subplot(1, M, m);
                imagesc(delay_axis * 1e3, time_axis, ...
                    20 * log10(squeeze(abs(u(:, m, :))).'), [-30, 0]);
                hold on;
                slide = -excess(m) * T;
                plot((tau0 - excess(m) * time_axis) * 1e3, time_axis, 'w--', ...
                    'LineWidth', 1);
                % The slide over the record, from the first arrival's start to
                % its end.
                plot([tau0, tau0 + slide] * 1e3, [T, T] * 0.92, 'w-', ...
                    'LineWidth', 2);
                plot(tau0 * 1e3, T * 0.92, 'w|', 'MarkerSize', 12, 'LineWidth', 2);
                plot((tau0 + slide) * 1e3, T * 0.92, 'w>', 'MarkerFaceColor', 'w');
                if slide < 0
                    set(findobj(gca, 'Marker', '>'), 'Marker', '<');
                end
                text((tau0 + slide / 2) * 1e3, T * 0.86, ...
                    sprintf('%+.2f ms', slide * 1e3), 'Color', 'w', ...
                    'HorizontalAlignment', 'center', 'FontWeight', 'bold');
                xlabel('Delay [ms]'); ylabel('Time [s]');
                title(sprintf('Hydrophone %d: f_{resamp} = 1/(1%+.1e)', m, excess(m)));
                colorbar;
            end
            sgtitle(sprintf(['Same CIR, different f_{resamp}: arrivals slide ', ...
                'by -(1/f_{resamp} - 1) T over T = %.2f s'], T));

            channel.f_resamp = [1, 1];
            testCase.verifyError(@() unpack(fs_out, 1, channel, buf, buf), ...
                'unpack:invalidFResamp');
        end

        function testInvalidInputs(testCase)
            % validate_inputs rejects each bad input with its own error
            % identifier, and lets a valid call through.

            L = 64; N_time = 50; M = 3;
            channel = struct();
            channel.h_hat = zeros(L, M, N_time);
            channel.h_hat(10, :, :) = 1;
            channel.version = 1.0;
            channel.params = struct('fs_delay', 8e3, 'fs_time', 100, 'fc', 12e3);
            fs = 100;

            testCase.verifyWarningFree(@() unpack(fs, 1:M, channel));
            testCase.verifyWarningFree(@() unpack(fs, [3, 1], channel, 0, 0.5));

            bad = channel; bad.version = 0.9;
            testCase.verifyError(@() unpack(fs, 1, bad), 'unpack:invalidVersion');

            for f = {0, -100, NaN, Inf, [100, 200], 1j}
                testCase.verifyError(@() unpack(f{1}, 1, channel), ...
                    'unpack:invalidFs', sprintf('fs = %s', mat2str(f{1})));
            end

            for index = {[], 0, M + 1, [1, 1], 1.5}
                testCase.verifyError(@() unpack(fs, index{1}, channel), ...
                    'unpack:invalidArrayIndex', ...
                    sprintf('array_index = %s', mat2str(index{1})));
            end

            testCase.verifyError(@() unpack(fs, 1, channel, 0.1), ...
                'unpack:invalidBuffer', 'one buffer only');
            for b = {-0.1, NaN, Inf, [0.1, 0.2]}
                testCase.verifyError(@() unpack(fs, 1, channel, b{1}, 0.1), ...
                    'unpack:invalidBuffer', sprintf('BUF_LEFT = %s', mat2str(b{1})));
                testCase.verifyError(@() unpack(fs, 1, channel, 0.1, b{1}), ...
                    'unpack:invalidBuffer', sprintf('BUF_RIGHT = %s', mat2str(b{1})));
            end

            for f = {[1, 1], [1, -1, 1], [1, NaN, 1], 0}
                bad = channel; bad.f_resamp = f{1};
                testCase.verifyError(@() unpack(fs, 1, bad), ...
                    'unpack:invalidFResamp', sprintf('f_resamp = %s', mat2str(f{1})));
            end
        end

        function testPhiHatTakesPrecedence(testCase)
            % A channel carrying both phi_hat and theta_hat unpacks as with
            % phi_hat alone.
            %
            % est saves its one phase trajectory as phi_hat with delay
            % tracking and as theta_hat without, so phi_hat is already the
            % whole phase and adding theta_hat would count the phase twice.
            % replay uses phi_hat alone too.

            fc = 15e3; fs_delay = 16e3; fs_time = 100; fs_out = 160;
            L = 128; N_time = 200;
            taps = zeros(L, 1);
            taps([21, 40, 63]) = [1; 0.5 - 0.2j; 0.25 + 0.1j];

            channel = struct();
            channel.h_hat = repmat(taps, [1, 1, N_time]);
            channel.version = 1.0;
            channel.params = struct('fs_delay', fs_delay, 'fs_time', fs_time, 'fc', fc);
            t = (0:ceil(N_time * fs_delay / fs_time) - 1) / fs_delay;
            channel.phi_hat = -2 * pi * fc * 1e-5 * t;
            u_phi = unpack(fs_out, 1, channel, 0.3, 0.3);
            channel.theta_hat = 0.7 * sin(2 * pi * 0.5 * t);
            u_both = unpack(fs_out, 1, channel, 0.3, 0.3);

            testCase.verifyEqual(u_both, u_phi, ['theta_hat changed the ', ...
                'unpacked channel although phi_hat is present.']);
        end

    end
end


%% Helpers

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


function samples = randsamples(population, num)
rand_index = randperm(length(population));
samples = population(rand_index(1:num));
end
