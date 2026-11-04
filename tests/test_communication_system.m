function tests = test_communication_system
% Focused checks for symbol recovery, real-noise calibration, and RF PSD units.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
projectRoot = fileparts(fileparts(mfilename('fullpath')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(projectRoot));
testCase.TestData.coefficients = load(fullfile(projectRoot, 'my_filter_coeffs_FIR.mat'), 'Num');
end

function setup(testCase)
previousRng = rng(42, 'twister');
testCase.addTeardown(@rng, previousRng);
end

function testConstellationMapping(testCase)
tx = Transmitter(4, 3.84e6, 1.92e9, 1, 'RRC');
tx.symbols = 0:3;
tx.do_mapping();
testCase.verifyEqual(tx.mapping.signal, [-1.5 -0.5 0.5 1.5]);
end

function testNoiselessRecovery(testCase)
for M = [2 4 8 16]
    [payload, received, firstIndex] = transmitFrame(testCase, M, Inf);
    testCase.assertGreaterThanOrEqual(firstIndex, 1);
    testCase.assertLessThanOrEqual(firstIndex + length(payload) - 1, length(received));
    testCase.verifyEqual(received(firstIndex:firstIndex+length(payload)-1), payload, ...
        sprintf('Noiseless %d-ASK payload must be recovered exactly.', M));
end
end

function testRealRfNoiseCalibration(testCase)
for snrDb = [6 12]
    [payload, received, firstIndex] = transmitFrame(testCase, 4, snrDb);
    testCase.assertGreaterThanOrEqual(firstIndex, 1);
    testCase.assertLessThanOrEqual(firstIndex + length(payload) - 1, length(received));
    ser = mean(payload ~= received(firstIndex:firstIndex+length(payload)-1));
    theory = 1.5 * qfunc(sqrt(6 * 10^(snrDb/10) / 15));
    % Allow finite-frame, synchronization, and blind-AGC differences.
    tolerance = 6 * sqrt(theory * (1-theory) / length(payload)) + 0.025;
    testCase.verifyEqual(ser, theory, 'AbsTol', tolerance);
end
end

function testPsdPowerAndFrequency(testCase)
oldVisibility = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
testCase.addTeardown(@set, groot, 'DefaultFigureVisible', oldVisibility);

Fs = 4.5e9;
fc = 1.92e9;
tx = Transmitter(4, 3.84e6, fc, 1, 'RRC');
tx.modulation.time = (0:8191) / Fs;
tx.modulation.signal = cos(2*pi*fc*tx.modulation.time);

for nfft = [4095 4096]
    tx.plot_modulator_psd(nfft);
    fig = gcf;
    testCase.addTeardown(@close, fig);
    ax = findobj(fig, 'Type', 'axes');
    curve = findobj(ax, 'Type', 'line');
    frequency = curve.XData * 1e6;
    density = 10.^(curve.YData/10);
    % A unit-amplitude real sinusoid has total mean-square power 1/2.
    testCase.verifyEqual(sum(density) * Fs/nfft, 0.5, 'RelTol', 0.01);
    [~, peak] = max(density);
    testCase.verifyEqual(abs(frequency(peak)), fc, 'AbsTol', Fs/nfft);
    testCase.verifyEqual(length(frequency), nfft);
    testCase.verifyTrue(contains(ax.YLabel.String, 'dB/Hz'));
end
end

function [payload, received, firstIndex] = transmitFrame(testCase, M, snrDb)
Rs = 3.84e6;
fc = 1.92e9;
Fs = 4.5e9;
Ns = 64;
barker = repmat([1 1 1 1 1 -1 -1 1 1 -1 1 -1 1], 1, 4);
preamble = (barker + 1) / 2 * (M - 1);
payload = randi([0 M-1], 1, 2000);

tx = Transmitter(M, Rs, fc, 1, 'RRC');
tx.symbols = [preamble payload];
tx.do_mapping();
tx.do_upsampling(Ns);
tx.do_pulse_shaping();
tx.do_dac(Fs);
tx.do_modulation();
signal = tx.modulation.signal;
if isfinite(snrDb)
    signal = awgn(signal, snrDb - 10*log10(Fs/(2*Rs)), 'measured');
end

rx = Receiver(M, Rs, fc, 1, Ns, testCase.TestData.coefficients.Num);
rx.do_demodulation(signal, tx.modulation.time);
rx.do_lowpass_filter();
rx.do_downsampling(tx.modulation.time);
rx.do_matched_filtering();
rx.do_detection();
received = rx.detection.symbols;
firstIndex = finddelay(preamble, received) + length(preamble) + 1;
end
