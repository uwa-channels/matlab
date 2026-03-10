classdef testNoise < matlab.unittest.TestCase
  % TESTNOISE Unit tests for noisegen.m
  %
  % Tests three noise generation modes:
  %   Option 1: Independent pink Gaussian noise (17 dB/decade).
  %             Invoked with nargin == 2.
  %   Option 2: Colored, spatially-correlated Gaussian noise
  %             (alpha = 2, beta mixing coefficients).
  %   Option 3: Impulsive (alpha-stable) noise
  %             (alpha < 2, beta mixing coefficients).
  %
  % Options 2 and 3 share a unified noise struct with fields:
  %   Fs, R, alpha, beta, fc, rms_power, version.
  % The stability index alpha determines the distribution:
  %   alpha == 2  => Gaussian (stabrnd reduces to Box-Muller).
  %   alpha <  2  => Symmetric alpha-stable (impulsive).
  % Bandpass filtering (fc +/- R/2*1.1) and rms_power scaling
  % are always applied when the noise struct is provided.
  %
  % The synthetic noise struct is modeled on real experimental
  % data (12-element ULA, 65 mixing taps, Fs=39062.5 Hz,
  % fc=13 kHz, R~4.9 kHz).  The beta mixing matrix is
  % constructed as a bandpass FIR with exponential spatial decay.
  %
  % Other m-files required: noisegen.m
  % Subfunctions: None
  % Toolbox required: Signal Processing Toolbox (for pwelch, bandpass,
  %                   fir1).
  % MAT-files required: None.
  %
  % See also: noisegen.m, replay.m
  %
  % Author: Zhengnan Li, Claude (Opus 4.6)
  % Email : uwa-channels@ofdm.link
  %
  % License: MIT
  %
  % Revision history:
  %   - Feb. 27, 2026: Initial release.
  %   - Mar.  9, 2026: Unified noise struct. Synthetic beta modeled
  %                     on real experimental data.

  properties (Constant)
    fs = 48000;
    N = 500000;        % For spectral / correlation / distribution tests
    N_SHORT = 100000;  % For size, resampling, subset, bandpass tests
  end

  methods (TestMethodSetup)
    function setRandomSeed(~)
      rng(1994);
    end
  end

  methods (Test)

    %% ============================================================
    %  Option 1: Pink Gaussian noise (nargin == 2)
    % =============================================================

    function testOption1_size(testCase)
      % Output size matches input size
      sizes = {[100000, 1], [100000, 4], [50000, 8]};
      for k = 1:length(sizes)
        w = noisegen(sizes{k}, testCase.fs);
        testCase.verifySize(w, sizes{k});
        testCase.verifyTrue(all(isfinite(w), 'all'));
      end
    end

    function testOption1_spectral_slope(testCase)
      % Verify ~17 dB/decade slope of pink noise PSD
      w = noisegen([testCase.N, 1], testCase.fs);
      [pxx, f] = pwelch(w, hann(8192), 4096, 8192, testCase.fs);

      % Fit slope between 100 Hz and 10 kHz
      mask = (f >= 100) & (f <= 10000);
      p = polyfit(log10(f(mask)), 10*log10(pxx(mask)), 1);

      % Plot
      figure('Name', 'Option1_PSD');
      subplot(121);
      plot(log10(f(mask)), 10*log10(pxx(mask)), 'DisplayName', 'Estimated'); hold on;
      plot(log10(f(mask)), polyval(p, log10(f(mask))), 'r--', 'DisplayName', ...
        sprintf('Fit: %.1f dB/dec', p(1)));
      psd_true = -17 * log10(f(mask)/1e3);
      psd_true = psd_true - mean(psd_true) + mean(10*log10(pxx(mask)));
      plot(log10(f(mask)), psd_true, 'k--', 'DisplayName', 'True (-17 dB/dec)');
      xlabel('log_{10}(f)'); ylabel('PSD [dB]');
      title('Pink noise PSD');
      legend; grid on;

      subplot(122);
      histogram(w, 100, 'Normalization', 'pdf');
      xlabel('Amplitude'); ylabel('PDF');
      title('Amplitude distribution');

      % Slope should be approximately -17 dB/decade
      testCase.verifyEqual(p(1), -17, 'RelTol', 0.15, ...
        sprintf('Pink noise slope %.1f dB/decade, expected ~-17', p(1)));
    end

    function testOption1_spatial_independence(testCase)
      % Channels should be independent
      w = noisegen([testCase.N, 4], testCase.fs);
      C = corrcoef(w);
      off_diag = C - eye(4);

      % Cross-correlation should be near zero
      testCase.verifyLessThan(max(abs(off_diag(:))), 0.1, ...
        'Channels are not independent');
    end

    %% ============================================================
    %  Option 2: Gaussian noise via beta mixing (alpha = 2)
    % =============================================================

    function testOption2_size_and_finite(testCase)
      noise = make_noise_struct(2);
      M = size(noise.beta, 1);
      input_size = [testCase.N_SHORT, M];
      w = noisegen(input_size, testCase.fs, 1:M, noise);
      testCase.verifySize(w, input_size);
      testCase.verifyTrue(all(isfinite(w), 'all'));
    end

    function testOption2_spatial_correlation(testCase)
      % Verify that the sample correlation matches the theoretical
      % prediction from the mixing equation:
      %   R_{ii'}(0) = sum_j sum_k beta(i,j,k) * beta(i',j,k)
      % i.e. R(0) = sum_k B_k * B_k^T, then normalize to get C.
      %
      % Bandpass is removed here to isolate the mixing test;
      % the bandpass filter modifies the cross-spectral density
      % non-uniformly across channel pairs, which breaks the
      % simple R = sum B_k B_k^T formula.  Bandpass correctness
      % Bandpass filtering preserves the correlation structure
      % only when all beta(i,j,:) profiles are proportional to
      % the same base FIR.  Setting perturb=0 ensures this.
      % Bandpass correctness is verified in testOption2_bandpass.
      noise = make_noise_struct(2, 0);
      M = size(noise.beta, 1);
      K = size(noise.beta, 3);
      input_size = [testCase.N, M];
      w = noisegen(input_size, testCase.fs, 1:M, noise);

      % --- Theoretical correlation from beta ---
      R_theory = zeros(M);
      for k = 1:K
        Bk = noise.beta(:, :, k);
        R_theory = R_theory + Bk * Bk.';
      end
      d = sqrt(diag(R_theory));
      C_theory = R_theory ./ (d * d.');

      % --- Sample correlation ---
      C_sample = corrcoef(w);

      % --- Plot ---
      figure('Name', 'Option2_correlation');

      subplot(221); hold on;
      for m = 1:M
        [pxx, f] = pwelch(w(:, m), hann(8192), 4096, 8192, testCase.fs);
        plot(f/1e3, 10*log10(pxx), 'DisplayName', sprintf('Ch %d', m));
      end
      xlabel('Frequency [kHz]'); ylabel('PSD [dB]');
      title('Estimated PSD'); legend('Location', 'best'); grid on;

      subplot(222);
      imagesc(C_theory); colorbar; axis square;
      title('Theoretical C'); xlabel('Channel'); ylabel('Channel');

      subplot(223);
      imagesc(C_sample); colorbar; axis square;
      title('Sample C'); xlabel('Channel'); ylabel('Channel');

      subplot(224);
      imagesc(C_sample - C_theory); colorbar; axis square;
      title('Error (sample - theory)'); xlabel('Channel'); ylabel('Channel');

      err = max(abs(C_sample(:) - C_theory(:)));
      sgtitle(sprintf('Gaussian mixing: max |error| = %.4f', err));

      % Max absolute error should be small
      testCase.verifyLessThan(err, 0.05, ...
        sprintf('Correlation mismatch: max |error| = %.4f', err));
    end

    function testOption2_identity_beta_independence(testCase)
      % Identity mixing should yield independent channels
      M = 12;
      noise = make_noise_struct(2);
      noise.beta = repmat(eye(M), [1, 1, 65]);

      input_size = [testCase.N_SHORT, M];
      w = noisegen(input_size, testCase.fs, 1:M, noise);
      C = corrcoef(w);
      off_diag = C - eye(M);

      testCase.verifyLessThan(max(abs(off_diag(:))), 0.1, ...
        'Diagonal-only beta should produce independent channels');
    end

    function testOption2_resampling(testCase)
      % Correct output length when fs != noise.Fs
      noise = make_noise_struct(2);
      M = size(noise.beta, 1);
      input_size = [testCase.N_SHORT, M];
      w = noisegen(input_size, testCase.fs, 1:M, noise);
      testCase.verifySize(w, input_size);
    end

    function testOption2_array_index_subset(testCase)
      % Using a subset of array indices
      noise = make_noise_struct(2);
      array_index = [2, 6, 10];
      input_size = [100000, length(array_index)];
      w = noisegen(input_size, testCase.fs, array_index, noise);
      testCase.verifySize(w, input_size);
      testCase.verifyTrue(all(isfinite(w), 'all'));
    end

    function testOption2_rms_scaling(testCase)
      % Verify per-channel rms_power scaling.
      % Use identity beta to avoid cross-channel leakage.
      noise = make_noise_struct(2);
      M = size(noise.beta, 1);
      noise.beta = repmat(eye(M), [1, 1, 65]);
      noise.rms_power = ones(M, 1);
      noise.rms_power(1) = 3;
      noise.rms_power(2) = 1;

      input_size = [500000, 2];
      w = noisegen(input_size, testCase.fs, [1, 2], noise);
      rms_ratio = rms(w(:, 1)) / rms(w(:, 2));

      testCase.verifyEqual(rms_ratio, 3, 'RelTol', 0.15, ...
        sprintf('RMS ratio %.2f, expected ~3', rms_ratio));
    end

    function testOption2_bandpass(testCase)
      % Verify bandpass filtering
      noise = make_noise_struct(2);
      fc = noise.fc;
      R = noise.R;

      input_size = [testCase.N_SHORT, 2];
      w = noisegen(input_size, testCase.fs, [1, 2], noise);

      [pxx, f] = pwelch(w(:, 1), hann(8192), 4096, 8192, testCase.fs);
      pxx_dB = 10*log10(pxx);

      in_band = (f >= fc - R/2) & (f <= fc + R/2);
      out_low = (f >= 100) & (f <= fc - R);
      out_high = (f >= fc + R) & (f <= testCase.fs/2 - 100);

      psd_in = mean(pxx_dB(in_band));
      psd_out = mean(pxx_dB(out_low | out_high));

      % Plot
      figure('Name', 'Option2_bandpass');
      plot(f/1e3, pxx_dB); hold on;
      xline((fc - R/2)/1e3, 'r--'); xline((fc + R/2)/1e3, 'r--');
      xlabel('Frequency [kHz]'); ylabel('PSD [dB]');
      title(sprintf('Gaussian bandpass: fc=%.0f kHz, R=%.1f kHz', fc/1e3, R/1e3));
      grid on;

      testCase.verifyGreaterThan(psd_in - psd_out, 20, ...
        sprintf('Bandpass rejection %.1f dB, expected > 20 dB', ...
        psd_in - psd_out));
    end

    function testOption2_gaussianity(testCase)
      % alpha=2 should produce Gaussian output (kurtosis ~ 3)
      noise = make_noise_struct(2);
      input_size = [testCase.N, 1];
      w = noisegen(input_size, testCase.fs, 1, noise);

      k = kurtosis(w);

      figure('Name', 'Option2_gaussianity');
      histogram(w, 200, 'Normalization', 'pdf'); hold on;
      x = linspace(min(w), max(w), 500);
      plot(x, normpdf(x, 0, std(w)), 'r', 'LineWidth', 1.5);
      xlabel('Amplitude'); ylabel('PDF');
      title(sprintf('Gaussianity check (kurtosis = %.2f)', k));

      testCase.verifyEqual(k, 3, 'RelTol', 0.15, ...
        sprintf('Kurtosis %.2f, expected ~3 for Gaussian', k));
    end

    %% ============================================================
    %  Option 3: Impulsive (alpha-stable) noise (alpha < 2)
    % =============================================================

    function testOption3_size_and_finite(testCase)
      noise = make_noise_struct(1.7);
      M = size(noise.beta, 1);
      input_size = [testCase.N_SHORT, M];
      w = noisegen(input_size, testCase.fs, 1:M, noise);
      testCase.verifySize(w, input_size);
      testCase.verifyTrue(all(isfinite(w), 'all'));
    end

    function testOption3_various_alpha(testCase)
      % Test a range of stability indices
      alphas = [1.2, 1.5, 1.7, 1.9];
      input_size = [100000, 2];
      for k = 1:length(alphas)
        noise = make_noise_struct(alphas(k));
        w = noisegen(input_size, testCase.fs, [1, 2], noise);
        testCase.verifySize(w, input_size);
        testCase.verifyTrue(all(isfinite(w), 'all'), ...
          sprintf('Non-finite values for alpha = %.1f', alphas(k)));
      end
    end

    function testOption3_heavier_tail(testCase)
      % Lower alpha should produce heavier tails (higher kurtosis)
      input_size = [500000, 1];
      noise_heavy = make_noise_struct(1.2);
      noise_light = make_noise_struct(1.9);

      rng(1994);
      w_heavy = noisegen(input_size, testCase.fs, 1, noise_heavy);
      rng(1994);
      w_light = noisegen(input_size, testCase.fs, 1, noise_light);

      % Plot
      figure('Name', 'Option3_tails');
      edges = linspace(-10, 10, 200);
      histogram(w_light, edges, 'Normalization', 'pdf', ...
        'DisplayName', '\alpha=1.9'); hold on;
      histogram(w_heavy, edges, 'Normalization', 'pdf', ...
        'DisplayName', '\alpha=1.2');
      xlabel('Amplitude'); ylabel('PDF');
      title('Tail comparison'); legend;
      set(gca, 'YScale', 'log');

      k_heavy = kurtosis(w_heavy);
      k_light = kurtosis(w_light);
      testCase.verifyGreaterThan(k_heavy, k_light, ...
        sprintf('alpha=1.2 kurtosis (%.1f) should exceed alpha=1.9 (%.1f)', ...
        k_heavy, k_light));
    end

    function testOption3_rms_scaling(testCase)
      % Verify rms_power scaling.
      % Use identity beta to avoid cross-channel leakage.
      noise = make_noise_struct(1.9);
      M = size(noise.beta, 1);
      noise.beta = repmat(eye(M), [1, 1, 65]);
      noise.rms_power = ones(M, 1);
      noise.rms_power(1) = 2;
      noise.rms_power(2) = 0.5;

      input_size = [500000, 2];
      w = noisegen(input_size, testCase.fs, [1, 2], noise);
      rms_ratio = rms(w(:, 1)) / rms(w(:, 2));

      testCase.verifyEqual(rms_ratio, 4, 'RelTol', 0.5, ...
        sprintf('RMS ratio %.2f, expected ~4', rms_ratio));
    end

    function testOption3_resampling(testCase)
      % Correct output when noise.Fs != fs
      noise = make_noise_struct(1.7);
      input_size = [100000, 2];
      w = noisegen(input_size, testCase.fs, [1, 2], noise);
      testCase.verifySize(w, input_size);
      testCase.verifyTrue(all(isfinite(w), 'all'));
    end

    function testOption3_bandpass(testCase)
      % Verify bandpass filtering for impulsive noise
      noise = make_noise_struct(1.7);
      fc = noise.fc;
      R = noise.R;

      input_size = [testCase.N_SHORT, 2];
      w = noisegen(input_size, testCase.fs, [1, 2], noise);

      [pxx, f] = pwelch(w(:, 1), hann(8192), 4096, 8192, testCase.fs);
      pxx_dB = 10*log10(pxx);

      in_band = (f >= fc - R/2) & (f <= fc + R/2);
      out_low = (f >= 100) & (f <= fc - R);
      out_high = (f >= fc + R) & (f <= testCase.fs/2 - 100);

      psd_in = mean(pxx_dB(in_band));
      psd_out = mean(pxx_dB(out_low | out_high));

      % Plot
      figure('Name', 'Option3_bandpass');
      plot(f/1e3, pxx_dB); hold on;
      xline((fc - R/2)/1e3, 'r--'); xline((fc + R/2)/1e3, 'r--');
      xlabel('Frequency [kHz]'); ylabel('PSD [dB]');
      title(sprintf('Impulsive bandpass: fc=%.0f kHz, R=%.1f kHz', ...
        fc/1e3, R/1e3));
      grid on;

      testCase.verifyGreaterThan(psd_in - psd_out, 20, ...
        sprintf('Bandpass rejection %.1f dB, expected > 20 dB', ...
        psd_in - psd_out));
    end

    function testOption3_alpha2_matches_gaussian(testCase)
      % alpha = 2 should look Gaussian-like (kurtosis ~ 3),
      % while alpha = 1.5 should have heavier tails
      input_size = [500000, 1];

      rng(42);
      noise_g = make_noise_struct(2);
      w_g = noisegen(input_size, testCase.fs, 1, noise_g);

      rng(42);
      noise_i = make_noise_struct(1.5);
      w_i = noisegen(input_size, testCase.fs, 1, noise_i);

      k_g = kurtosis(w_g);
      k_i = kurtosis(w_i);

      testCase.verifyGreaterThan(k_i, k_g, ...
        sprintf('alpha=1.5 kurtosis (%.1f) should exceed alpha=2 (%.1f)', ...
        k_i, k_g));
    end

    function testOption3_spatial_correlation(testCase)
      % For alpha < 2, Pearson correlation overestimates dependence
      % because heavy-tailed outliers inflate the sample covariance.
      % Spearman (rank) correlation is robust to this and tracks the
      % Gaussian-predicted structure more closely.
      %
      % With zero perturbation, all beta profiles are proportional
      % to the same base FIR, so the bandpass does not alter the
      % correlation structure.
      noise = make_noise_struct(1.7, 0);
      M = size(noise.beta, 1);
      K = size(noise.beta, 3);
      input_size = [testCase.N, M];
      w = noisegen(input_size, testCase.fs, 1:M, noise);

      C_pearson = corrcoef(w);
      % Spearman = Pearson on ranks (no Statistics Toolbox needed)
      w_ranked = zeros(size(w));
      for m = 1:size(w, 2)
        [~, idx] = sort(w(:, m));
        w_ranked(idx, m) = (1:size(w, 1)).';
      end
      C_spearman = corrcoef(w_ranked);

      % Gaussian-predicted correlation (reference)
      R_theory = zeros(M);
      for k = 1:K
        Bk = noise.beta(:, :, k);
        R_theory = R_theory + Bk * Bk.';
      end
      d = sqrt(diag(R_theory));
      C_theory = R_theory ./ (d * d.');

      % Plot
      figure('Name', 'Option3_correlation');

      subplot(221);
      imagesc(C_theory); colorbar; axis square;
      title('Gaussian-predicted C'); xlabel('Channel'); ylabel('Channel');

      subplot(222);
      imagesc(C_pearson); colorbar; axis square;
      title(sprintf('Pearson C (\\alpha=%.1f)', noise.alpha));
      xlabel('Channel'); ylabel('Channel');

      subplot(223);
      imagesc(C_spearman); colorbar; axis square;
      title(sprintf('Spearman C (\\alpha=%.1f)', noise.alpha));
      xlabel('Channel'); ylabel('Channel');

      subplot(224);
      offsets = 1:M-1;
      c_theory_d = arrayfun(@(d) mean(diag(C_theory, d)), offsets);
      c_pearson_d = arrayfun(@(d) mean(diag(C_pearson, d)), offsets);
      c_spearman_d = arrayfun(@(d) mean(diag(C_spearman, d)), offsets);
      plot(offsets, c_theory_d, 'ko-', 'DisplayName', 'Gaussian theory'); hold on;
      plot(offsets, c_pearson_d, 'rs-', 'DisplayName', 'Pearson');
      plot(offsets, c_spearman_d, 'b^-', 'DisplayName', 'Spearman');
      xlabel('Element separation |i-j|'); ylabel('Mean correlation');
      title('Correlation decay'); legend; grid on;

      sgtitle(sprintf('Impulsive (\\alpha=%.1f) spatial correlation', noise.alpha));

      % Structural checks (Spearman is the robust metric)
      testCase.verifyGreaterThan(abs(C_spearman(1,2)), abs(C_spearman(1, M)), ...
        'Adjacent channels should be more correlated than distant ones');

      c_vs_d = arrayfun(@(d) mean(abs(diag(C_spearman, d))), 1:M-1);
      testCase.verifyTrue(c_vs_d(1) > c_vs_d(end), ...
        'Correlation should decay with element separation');
    end

  end
end


%% ====================================================================
%  Helper functions
% =====================================================================

function noise = make_noise_struct(alpha, perturb)
% Construct a synthetic noise struct that mimics real experimental
% data from a 12-element ULA.
%
% The beta mixing matrix is built as follows:
%   1. Design a bandpass FIR (K=65 taps) at the channel center
%      frequency.  This captures the spectral shaping observed
%      in real measured mixing coefficients.
%   2. Scale each (i,j) pair by 0.5^|i-j| to model the spatial
%      correlation decay of a linear array.
%   3. Add random perturbations scaled by PERTURB (default 0.05),
%      as observed in real data.
%
% Parameters are modeled on a 12-element array with:
%   Fs = 39062.5 Hz, fc = 13 kHz, R ~ 4.9 kHz.

if nargin < 2
  perturb = 0.05;
end

M = 12;           % Number of array elements
K = 65;           % Number of mixing taps
Fs = 39062.5;     % Noise sampling rate [Hz]
fc = 13000;       % Center frequency [Hz]
R = 4882.8125;    % Bandwidth [Hz]

% --- Bandpass FIR as base tap profile ---
% Normalized cutoff frequencies for fir1
f_lo = (fc - R/2) / (Fs/2);
f_hi = (fc + R/2) / (Fs/2);
base_profile = fir1(K-1, [f_lo, f_hi]);
base_profile = base_profile(:);

% --- Build beta with spatial decay ---
rng_state = rng;
rng(2024);  % Deterministic for reproducibility

decay_rate = 0.5;

beta = zeros(M, M, K);
for i = 1:M
  for j = 1:M
    d = abs(i - j);
    scale = decay_rate^d;
    perturbation = perturb * scale * randn(K, 1);
    beta(i, j, :) = scale * base_profile + perturbation;
  end
end

rng(rng_state);  % Restore RNG state

% --- Per-channel RMS power ---
% Mild variation across elements (ratio max/min ~ 1.3)
rms_power = 3.2e-4 * (1 + 0.1 * linspace(-1, 1, M).');

noise = struct();
noise.Fs = Fs;
noise.R = R;
noise.alpha = alpha;
noise.beta = beta;
noise.fc = fc;
noise.rms_power = rms_power;
noise.version = 1.0;
end


function k = kurtosis(x)
x = x(:);
mu = mean(x);
m2 = mean((x - mu).^2);
m4 = mean((x - mu).^4);
k = m4 / m2^2;
end

function y = normpdf(x, mu, sigma)
y = 1 / (sigma * sqrt(2*pi)) * exp(-0.5 * ((x - mu) / sigma).^2);
end
