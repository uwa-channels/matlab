classdef testNoise < matlab.unittest.TestCase
  % TESTNOISE Unit tests for noisegen.m
  %
  % Tests three noise generation modes:
  %   Option 1: Independent pink Gaussian noise (17 dB/decade).
  %   Option 2: Colored, spatially-correlated Gaussian noise.
  %   Option 3: Impulsive (alpha-stable) noise.
  %
  % Verifies: output size, spectral shape, spatial correlation,
  % resampling, array indexing, bandpass filtering, and power scaling.
  %
  % Other m-files required: noisegen.m
  % Subfunctions: None
  % Toolbox required: Signal Processing Toolbox (for pwelch, bandpass).
  % MAT-files required: None.
  %
  % See also: noisegen.m, replay.m
  %
  % Author: Zhengnan Li
  % Email : uwa-channels@ofdm.link
  %
  % License: MIT
  %
  % Revision history:
  %   - Feb. 27, 2026: Initial release.

  properties (Constant)
    fs = 48000;
    N = 1000000;  % Long enough for spectral estimation
  end

  methods (TestMethodSetup)
    function setRandomSeed(~)
      rng(1994);
    end
  end

  methods (Test)

    %% ============================================================
    %  Option 1: Pink Gaussian noise
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
    %  Option 2: Colored spatially-correlated Gaussian noise
    % =============================================================

    function testOption2_size_and_finite(testCase)
      input_size = [testCase.N, 4];
      array_index = 1:4;
      noise = make_gaussian_noise(4, testCase.fs);
      w = noisegen(input_size, testCase.fs, array_index, noise);
      testCase.verifySize(w, input_size);
      testCase.verifyTrue(all(isfinite(w), 'all'));
    end

    function testOption2_spatial_correlation(testCase)
      % Verify that spatial correlation structure is preserved
      M = 4;
      % Create a non-trivial correlation matrix
      A = randn(M); A = A * A';
      sigma = A / max(abs(A(:)));  % normalize
      sigma = sigma + M * eye(M);  % ensure PD

      noise = struct();
      noise.Fs = testCase.fs;
      noise.sigma = sigma;
      noise.h = zeros(100, M); noise.h(50, :) = 1;  % delta filter
      noise.version = 1.0;

      input_size = [testCase.N, M];
      w = noisegen(input_size, testCase.fs, 1:M, noise);

      % Estimate sample correlation
      C_sample = corrcoef(w);
      C_truth = diag(1./sqrt(diag(sigma))) * sigma * diag(1./sqrt(diag(sigma)));

      % Should roughly match (within statistical tolerance)
      err = max(abs(C_sample(:) - C_truth(:)));
      testCase.verifyLessThan(err, 0.15, ...
        sprintf('Correlation mismatch: max error = %.3f', err));
    end

    function testOption2_resampling(testCase)
      % Verify correct output length when fs != noise.Fs
      noise_fs = 96000;
      input_size = [testCase.N, 2];
      noise = make_gaussian_noise(2, noise_fs);
      w = noisegen(input_size, testCase.fs, [1, 2], noise);
      testCase.verifySize(w, input_size);
    end

    function testOption2_array_index_subset(testCase)
      % Using a subset of array indices
      M_total = 6;
      noise = make_gaussian_noise(M_total, testCase.fs);
      array_index = [2, 4, 6];
      input_size = [100000, length(array_index)];
      w = noisegen(input_size, testCase.fs, array_index, noise);
      testCase.verifySize(w, input_size);
      testCase.verifyTrue(all(isfinite(w), 'all'));
    end

    %% ============================================================
    %  Option 3: Impulsive (alpha-stable) noise
    % =============================================================

    function testOption3_size_and_finite(testCase)
      input_size = [testCase.N, 3];
      array_index = 1:3;
      noise = make_impulsive_noise(1.7, 3, testCase.fs);
      w = noisegen(input_size, testCase.fs, array_index, noise);
      testCase.verifySize(w, input_size);
      testCase.verifyTrue(all(isfinite(w), 'all'));
    end

    function testOption3_various_alpha(testCase)
      % Test a range of stability indices
      alphas = [1.2, 1.5, 1.7, 1.9];
      input_size = [100000, 2];
      for k = 1:length(alphas)
        noise = make_impulsive_noise(alphas(k), 2, testCase.fs);
        w = noisegen(input_size, testCase.fs, [1, 2], noise);
        testCase.verifySize(w, input_size);
        testCase.verifyTrue(all(isfinite(w), 'all'), ...
          sprintf('Non-finite values for alpha = %.1f', alphas(k)));
      end
    end

    function testOption3_heavier_tail(testCase)
      % Lower alpha should produce heavier tails (more outliers)
      input_size = [500000, 1];
      noise_heavy = make_impulsive_noise(1.2, 1, testCase.fs);
      noise_light = make_impulsive_noise(1.9, 1, testCase.fs);

      rng(1994);
      w_heavy = noisegen(input_size, testCase.fs, 1, noise_heavy);
      rng(1994);
      w_light = noisegen(input_size, testCase.fs, 1, noise_light);

      % Kurtosis should be higher for heavier tails
      k_heavy = kurtosis(w_heavy);
      k_light = kurtosis(w_light);
      testCase.verifyGreaterThan(k_heavy, k_light, ...
        sprintf('alpha=1.2 kurtosis (%.1f) should exceed alpha=1.9 (%.1f)', ...
        k_heavy, k_light));
    end

    function testOption3_rms_scaling(testCase)
      % Verify rms_power scaling
      input_size = [200000, 2];
      noise = make_impulsive_noise(1.7, 2, testCase.fs);
      noise.rms_power = [2; 0.5];

      w = noisegen(input_size, testCase.fs, [1, 2], noise);
      rms_ratio = rms(w(:, 1)) / rms(w(:, 2));

      % Should be approximately 2/0.5 = 4
      testCase.verifyEqual(rms_ratio, 4, 'RelTol', 0.3, ...
        sprintf('RMS ratio %.2f, expected ~4', rms_ratio));
    end

    function testOption3_bandpass(testCase)
      % Verify bandpass filtering when fc and R are specified
      fc = 12000;
      R = 4000;
      input_size = [testCase.N, 1];
      noise = make_impulsive_noise(1.7, 1, testCase.fs);
      noise.fc = fc;
      noise.R = R;

      w = noisegen(input_size, testCase.fs, 1, noise);
      [pxx, f] = pwelch(w, hann(8192), 4096, 8192, testCase.fs);

      % Power should be concentrated in [fc-R/2, fc+R/2]
      in_band = (f >= fc - R/2) & (f <= fc + R/2);
      out_band = (f < fc - R) | (f > fc + R);
      out_band(f == 0) = false;

      pwr_in = mean(pxx(in_band));
      pwr_out = mean(pxx(out_band));
      rejection_dB = 10*log10(pwr_in / pwr_out);

      testCase.verifyGreaterThan(rejection_dB, 15, ...
        sprintf('Bandpass rejection %.1f dB, expected > 15 dB', rejection_dB));
    end

    function testOption3_resampling(testCase)
      % Verify correct output when noise.Fs != fs
      noise_fs = 96000;
      input_size = [100000, 2];
      noise = make_impulsive_noise(1.7, 2, noise_fs);
      w = noisegen(input_size, testCase.fs, [1, 2], noise);
      testCase.verifySize(w, input_size);
      testCase.verifyTrue(all(isfinite(w), 'all'));
    end

  end
end


%% ====================================================================
%  Helper functions
% =====================================================================

function noise = make_gaussian_noise(M, Fs)
noise = struct();
noise.Fs = Fs;
noise.sigma = eye(M);
noise.h = randn(100, M);
noise.version = 1.0;
end

function noise = make_impulsive_noise(alpha, M, Fs)
noise = struct();
noise.alpha = alpha;
noise.Fs = Fs;
noise.beta = repmat(eye(M), [1, 1, 65]);
noise.version = 1.0;
end

function k = kurtosis(x)
x = x(:);
mu = mean(x);
m2 = mean((x - mu).^2);
m4 = mean((x - mu).^4);
k = m4 / m2^2;
end
