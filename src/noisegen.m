function w = noisegen(input_size, fs, varargin)
% NOISEGEN Generate underwater acoustic noise.
%
% W = NOISEGEN(INPUT_SIZE, FS) generates independent pink Gaussian
% noise with 17 dB / decade power spectrum, W. INPUT_SIZE specifies the
% size of the input signal, and FS is the sampling rate of the input signal.
% The noise is assumed to be sptailly independent.
%
% W = NOISEGEN(INPUT_SIZE, FS, ARRAY_INDEX, NOISE) generates acoustic
% noise. Based on the input of the NOISE struct, it can generate colored
% and spatially-correlated Gaussian noise or impulsive noise. NOISE is a
% struct containing noise statistics, and ARRAY_INDEX specifies the
% indices of the receiver elements through which the user would like to
% generate noise. These indices must match those used in the REPLAY function.
%
% Inputs:
%    input_size              - Input signal matrix size.
%    fs                      - Sampling rate of the input signal in Hz.
%    array_index             - Indices of the channels.
%    noise                   - Struct containing noise statistics.
%
% Outputs:
%    w                       - Generated noise.
%
% Example:
%    See example_replay.m
%
%
% Other m-files required: None
% Subfunctions: validate_inputs, noise_option_1, noise_option_2, noise_option_3
% Toolbox required: Signal Processing Toolbox (R) (resample function).
% MAT-files optional: noise matfile.
%
% See also: replay.m
%
% Authors: Zhengnan Li, Diego A. Cuji, and Mandar Chitre
%
% Emails: uwa-channels@ofdm.link, cujidutan.d@northeastern.edu, mandar@nus.edu.sg
% License: MIT
%
% Revision history:
%   - Apr. 1, 2025: initial release.
%
%


if nargin ~= 2
    array_index = varargin{1};
    noise = varargin{2};
    validate_inputs(input_size, noise, array_index);
end


if nargin == 2
    w = noise_option_1(input_size, fs);
elseif nargin == 4
    if ~isfield(noise, 'alpha')
        w = noise_option_2(input_size, fs, noise, array_index);
    else
        w = noise_option_3(input_size, fs, noise, array_index);
    end
else
    error('Wrong noise option.');
end
end

function validate_inputs(input_size, noise, array_index)
% Noise version checking
assert(noise.version >= 1.0, ...
    'The minimum version of the noise matrix is 1.0, and you have %.1f.', ...
    noise.version);

if isfield(noise, 'alpha')
    M = size(noise.beta, 1);
else
    M = size(noise.sigma, 1);
end

assert(length(unique(array_index)) == length(array_index));
assert(input_size(2) <= M, 'The number of signal channels, %d, should be less than %d.', ...
    input_size(2), M)
assert(max(array_index) <= M, ...
    'The largest receive channel array_index, %d, should be less than %d.', ...
    max(array_index), M);
end


function w = noise_option_1(input_size, fs)
% Generate textbook style noise: independent pink Gaussian noise (17dB/decade) across array elements.
if mod(input_size(1), 2) == 1
    N = input_size(1) + 1;
else
    N = input_size(1);
end
M = input_size(2);
white_noise = randn(N, M);
spectrum = fft(white_noise, [], 1);
freqs = (0:N / 2)' * fs / N;
alpha = 1.7;
half_N = floor(N/2);
filter_shape = ones(N, M);
filter_shape(1:half_N, :) = repmat(1./(freqs(2:half_N+1).^(alpha / 2)), 1, M);
filter_shape(half_N:end-1, :) = flip(filter_shape(1:half_N, :), 1);
filtered_spectrum = spectrum .* filter_shape;
w = real(ifft(filtered_spectrum));
end


function w = noise_option_2(input_size, fs, noise, array_index)
% Generate noise according to statistics collected during experiments.
[p, q] = rat(fs/noise.Fs);
signal_size = input_size;
signal_size(1) = ceil(signal_size(1)*q/p);
n = randn(signal_size(1), size(noise.sigma, 1)) * chol(noise.sigma);
w = zeros(signal_size);
for m = 1:signal_size(2)
    w(:, m) = conv(n(:, array_index(m)), noise.h(:, array_index(m)), 'same');
end
w = w - mean(w, 1);
w_resampled = zeros(ceil(size(w, 1)*p/q), size(w, 2));
for m = 1:size(w, 2)
    w_resampled(:, m) = resample(w(:, m), p, q);
end
w = w(1:input_size(1), :);
end


function w = noise_option_3(input_size, fs, noise, array_index)
% Generate impulsive noise.
alpha = noise.alpha;
beta = noise.beta;
Fs = noise.Fs;

[p, q] = rat(fs/Fs);
signal_size = input_size;
signal_size(1) = ceil(signal_size(1)*q/p);

K = signal_size(1);
N = size(beta, 1);

z = stabrnd(alpha, 0, 1, 0, K+size(beta, 3), N);
w = zeros(K, N);
for i = 1:N
    for j = 1:N
        for k = 1:size(beta, 3)
            w(:, i) = w(:, i) + beta(i, j, k) .* z(k:k+K-1, j);
        end
    end
end

w = w(:, array_index);
w_resampled = zeros(ceil(size(w, 1)*p/q), size(w, 2));
for m = 1:size(w, 2)
    w_resampled(:, m) = resample(w(:, m), p, q);
end
w = w(1:input_size(1), :);
end


function x = stabrnd(alpha, beta, c, delta, m, n)

% STABRND Stable Random Number Generator
% (McCulloch 12/18/96)
%
% x = stabrnd(alpha,beta,c,delta,m,n);
%
% alpha, beta, c and delta are the characteristic exponent,
%   symmetry paramter, scale parameter (dispersion^1/alpha) and
%   location parameter respectively.
%
% Returns m x n matrix of iid stable random numbers with
%   characteristic exponent alpha in [.1,2], skewness parameter
%   beta in [-1,1], scale c > 0, and location parameter delta.
%
% Based on the method of J.M. Chambers, C.L. Mallows and B.W.
%   Stuck, "A Method for Simulating Stable Random Variables,"
%   JASA 71 (1976): 340-4.
%
% Encoded in MATLAB by J. Huston McCulloch, Ohio State
%   University Econ. Dept. (mcculloch.2@osu.edu).  This 12/18/96
%   version uses 2*m*n calls to RAND, and does not rely on
%   the STATISTICS toolbox.

% The CMS method is applied in such a way that x will have the
%   log characteristic function
%        log E exp(ixt) = i*delta*t + psi(c*t),
%   where
%     psi(t) = -abs(t)^alpha*(1-i*beta*sign(t)*tan(pi*alpha/2))
%                              for alpha ~= 1,
%            = -abs(t)*(1+i*beta*(2/pi)*sign(t)*log(abs(t))),
%                              for alpha = 1.
%
% With this parameterization, the stable cdf S(x; alpha, beta,
%   c, delta) equals S((x-delta)/c; alpha, beta, 1, 0).  See my
%   "On the parametrization of the afocal stable distributions,"
%   _Bull. London Math. Soc._ 28 (1996): 651-55, for details.
%
% When alpha = 2, the distribution is Gaussian with mean delta
%   and variance 2*c^2, and beta has no effect.
%
% When alpha > 1, the mean is delta for all beta.  When alpha
%   <= 1, the mean is undefined.
%
% When beta = 0, the distribution is symmetrical and delta is
%   the median for all alpha.  When alpha = 1 and beta = 0, the
%   distribution is Cauchy (arctangent) with median delta.
%
% When the submitted alpha is > 2 or < .1, or beta is outside
%   [-1,1], an error message is generated and x is returned as a
%   matrix of NaNs.
%
% Alpha < .1 is not allowed here because of the non-negligible
%   probability of overflows.
%
% If you're only interested in the symmetric cases, you may just
%   set beta = 0 and skip the following considerations:
%
% When beta > 0 (< 0), the distribution is skewed to the right
%   (left).
%
% When alpha < 1, delta, as defined above, is the unique fractile
%   that is invariant under averaging of iid contributions.  I
%   call such a fractile a "focus of stability."  This, like the
%   mean, is a natural location parameter.
%
% When alpha = 1, either every fractile is a focus of stability,
%   as in the beta = 0 Cauchy case, or else there is no focus of
%   stability at all, as is the case for beta ~=0.  In the latter
%   cases, which I call "afocal," delta is just an arbitrary
%   fractile that has a simple relation to the c.f.
%
% When alpha > 1 and beta > 0, med(x) must lie very far below
%   the mean as alpha approaches 1 from above.  Furthermore, as
%   alpha approaches 1 from below, med(x) must lie very far above
%   the focus of stability when beta > 0.  If beta ~= 0, there
%   is therefore a discontinuity in the distribution as a function
%   of alpha as alpha passes 1, when delta is held constant.
%
% CMS, following an insight of Vladimir Zolotarev, remove this
%   discontinuity by subtracting
%          beta*c*tan(pi*alpha/2)
%   (equivalent to their -tan(alpha*phi0)) from x for alpha ~=1
%   in their program RSTAB, a.k.a. RNSTA in IMSL (formerly GGSTA).
%   The result is a random number whose distribution is a contin-
%   uous function of alpha, but whose location parameter (which I
%   call zeta) is a shifted version of delta that has no known
%   interpretation other than computational convenience.
%   The present program restores the more meaningful "delta"
%   parameterization by using the CMS (4.1), but with
%   beta*c*tan(pi*alpha/2) added back in (ie with their initial
%   tan(alpha*phi0) deleted).  RNSTA therefore gives different
%   results than the present program when beta ~= 0.  However,
%   the present beta is equivalent to the CMS beta' (BPRIME).
%
% Rather than using the CMS D2 and exp2 functions to compensate
%   for the ill-condition of the CMS (4.1) when alpha is very
%   near 1, the present program merely fudges these cases by
%   computing x from their (2.4) and adjusting for
%   beta*c*tan(pi*alpha/2) when alpha is within 1.e-8 of 1.
%   This should make no difference for simulation results with
%   samples of size less than approximately 10^8, and then
%   only when the desired alpha is within 1.e-8 of 1, but not
%   equal to 1.
%
% The frequently used Gaussian and symmetric cases are coded
%   separately so as to speed up execution.
%
% Additional references:
% V.M. Zolotarev, _One Dimensional Stable Laws_, Amer. Math.
%   Soc., 1986.
% G. Samorodnitsky and M.S. Taqqu, _Stable Non-Gaussian Random
%   Processes_, Chapman & Hill, 1994.
% A. Janicki and A. Weron, _Simulaton and Chaotic Behavior of
%   Alpha-Stable Stochastic Processes_, Dekker, 1994.
% J.H. McCulloch, "Financial Applications of Stable Distributons,"
%   _Handbook of Statistics_ Vol. 14, forthcoming early 1997.

% Errortraps:
if alpha < .1 | alpha > 2
    disp('Alpha must be in [.1,2] for function STABRND.')
    alpha
    x = NaN * zeros(m, n);
    return
end
if abs(beta) > 1
    disp('Beta must be in [-1,1] for function STABRND.')
    beta
    x = NaN * zeros(m, n);
    return
end

% Generate exponential w and uniform phi:
w = -log(rand(m, n));
phi = (rand(m, n) - .5) * pi;

% Gaussian case (Box-Muller):
if alpha == 2
    x = (2 * sqrt(w) .* sin(phi));
    x = delta + c * x;
    return
end

% Symmetrical cases:
if beta == 0
    if alpha == 1 % Cauchy case
        x = tan(phi);
    else
        x = ((cos((1 - alpha)*phi) ./ w).^(1 / alpha - 1) ...
            .* sin(alpha * phi) ./ cos(phi).^(1 / alpha));
    end

    % General cases:
else
    cosphi = cos(phi);
    if abs(alpha-1) > 1.e-8
        zeta = beta * tan(pi*alpha/2);
        aphi = alpha * phi;
        a1phi = (1 - alpha) * phi;
        x = ((sin(aphi) + zeta * cos(aphi)) ./ cosphi) ...
            .* ((cos(a1phi) + zeta * sin(a1phi)) ...
            ./ (w .* cosphi)).^((1 - alpha) / alpha);
    else
        bphi = (pi / 2) + beta * phi;
        x = (2 / pi) * (bphi .* tan(phi) - beta * log((pi / 2)*w ...
            .*cosphi./bphi));
        if alpha ~= 1
            x = x + beta * tan(pi * alpha/2);
        end
    end
end

% Finale:
x = delta + c * x;
return
% End of STABRND.M
end

% [EOF]
