classdef testNoise < matlab.unittest.TestCase

    methods (TestMethodSetup)
        function setRandomSeed(~)
            rng(1994);
        end
    end

    methods (Test)

        function testNoiseOption1(testCase)
            input = zeros(1000, 4);
            fs = 48000;

            w = noisegen(size(input), fs);

            testCase.verifySize(w, size(input));
            testCase.verifyTrue(all(isfinite(w), 'all'));
        end

        function testNoiseOption2(testCase)
            input = zeros(1000, 4);
            array_index = [1, 2, 3, 4];
            fs = 48000;
            noise.Fs = fs;
            noise.sigma = eye(4);
            noise.h = randn(100, 4);
            noise.version = 1.0;

            w = noisegen(size(input), fs, array_index, noise);

            testCase.verifySize(w, size(input));
            testCase.verifyTrue(all(isfinite(w), 'all'));
        end

        function testPowerNormalization1(testCase)
            input = zeros(1e6, 4);
            fs = 48000;

            w = noisegen(size(input), fs);
            p = mean(abs(w).^2, 1);

            testCase.verifyTrue(all(abs(p-1) < 1e-3));
        end

        function testPowerNormalization2(testCase)
            input = zeros(1e6, 4);
            array_index = [1, 2, 3, 4];
            fs = 48000;
            noise.Fs = fs;
            noise.sigma = zeros(4);
            noise.sigma(1:5:end) = rand(1, 4);
            noise.h = randn(100, 4);
            noise.h = noise.h ./ max(abs(noise.h), 1);
            noise.version = 1.0;

            w = noisegen(size(input), fs, array_index, noise);
            p = mean(abs(w).^2, 1);
            metric = p.' ./ diag(noise.sigma);

            testCase.verifyTrue(all(abs(metric-metric(1)) < 1e-2));
        end

        function testImpulsiveNoise1(testCase)
            input = zeros(1000, 3);
            array_index = [1, 2, 3];
            fs = 48000;
            noise.alpha = 1.7;
            noise.Fs = 48000;
            noise.beta = repmat(diag([1, 2, 3]), [1, 1, 65]);
            noise.version = 1.0;

            w = noisegen(size(input), fs, array_index, noise);

            testCase.verifySize(w, size(input));
            testCase.verifyTrue(all(isfinite(w), 'all'));
        end

    end
end
