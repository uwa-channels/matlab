classdef testNoise < matlab.unittest.TestCase

    methods (TestMethodSetup)
        function setRandomSeed(~)
            rng(1994);
        end
    end

    methods (Test)

        function testNoiseOption1(testCase)
            input = zeros(1000000, 4);
            fs = 48000;

            w = noisegen(size(input), fs);

            testCase.verifySize(w, size(input));
            testCase.verifyTrue(all(isfinite(w), 'all'));
        end

        function testNoiseOption2(testCase)
            input = zeros(1000000, 4);
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

        function testNoiseOption3(testCase)
            input = zeros(1000000, 3);
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
