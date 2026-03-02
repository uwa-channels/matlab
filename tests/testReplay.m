classdef testReplay < matlab.unittest.TestCase
  % TESTREPLAY Unit tests for replay.m
  %
  % Tests {static, mobile} x {theta_hat, phi_hat} across various
  % speed configurations.
  %
  %   theta_hat: phase tracking only. h_hat drifts with motion.
  %   phi_hat:   phase + delay tracking. h_hat is static.
  %
  % Other m-files required: replay.m
  % Subfunctions: place_taps, compensate_doppler, randsamples
  % Toolbox required: Signal Processing Toolbox (for resample, findpeaks).
  % MAT-files required: None.
  %
  % See also: replay.m, est.m
  %
  % Author: Zhengnan Li, Claude (Opus 4.6)
  % Email : uwa-channels@ofdm.link
  %
  % License: MIT
  %
  % Revision history:
  %   - Feb. 27, 2026: Initial release.

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
      if p.n_cycles > 0
        phi = phi + p.fc * p.v_amp * T_ch ...
          / (c * p.n_cycles) ...
          * (cos(omega_osc * t_delay) - 1);
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
      xlabel('Time [s]'); ylabel('Speed [m/s]');
      title('Instantaneous speed (from \phi)');
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

        n_compare = min(n_found, p.n_path);
        tol = 2e-4 * p.n_path;
        if n_found > 0
          criteria(m) = abs( ...
            sum(est_delays(idx_e(1:n_compare)) .* est_gains(idx_e(1:n_compare)).') - ...
            sum(path_delay_0(idx_t(1:n_compare), m) .* path_gain(idx_t(1:n_compare), m))) ...
            < tol;
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
