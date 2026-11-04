%% M-ASK communication link: transmitter, AWGN channel, and receiver
% Run from the repository root. See README.md for requirements and parameters.
% The local Live Editor report is not required to run this source script.

clear; close all; clc;
rng(42, 'twister');
coefficients = load('my_filter_coeffs_FIR.mat', 'Num');

%% 1. System configuration
M = 4;                 % Modulation order: 2, 4, 8, or 16
Rs = 3.84e6;           % Symbol rate (Baud)
Ns = 64;               % Baseband samples per symbol
Vd = 1;                % Distance between adjacent voltage levels
fc = 1920e6;           % Carrier frequency (Hz)
Fs_analog = 4500e6;    % RF sample rate (Hz)
payload_len = 2000;    % Increase for more precise SER estimates (uses more RAM)
target_snr_db = 10;    % Es/N0 in dB, using the real-noise PSD convention N0/2

fprintf('System: %d-ASK | Rs = %.2f MBaud | Fs = %.2f GHz | fc = %.2f GHz\n', ...
    M, Rs/1e6, Fs_analog/1e9, fc/1e9);

%% 2. Barker preamble and random payload
barker13 = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
seq_binary = repmat(barker13, 1, 4);
preamble = zeros(size(seq_binary));
preamble(seq_binary == 1) = M - 1;
preamble_len = length(preamble);
payload = randi([0 M-1], 1, payload_len);

fprintf('Frame: %d preamble symbols + %d payload symbols\n', ...
    preamble_len, payload_len);

%% 3. Transmitter and intermediate signal plots
tx = Transmitter(M, Rs, fc, Vd, 'RRC');
tx.symbols = [preamble, payload];
tx.do_mapping();
tx.do_upsampling(Ns);
tx.do_pulse_shaping();
tx.do_dac(Fs_analog);
tx.do_modulation();

s_t = tx.modulation.signal;
t_axis = tx.modulation.time;
tx.plot_constellation();
tx.plot_baseband_signals(70);
tx.plot_modulator_output(70 / Rs);
tx.plot_modulator_psd();

%% 4. Real-valued RF AWGN channel
% Es/N0 = Fs_analog/(2*Rs) * SNR_sample for real noise with PSD N0/2.
process_gain = 10 * log10(Fs_analog / (2 * Rs));
input_snr = target_snr_db - process_gain;
r_t = awgn(s_t, input_snr, 'measured');
fprintf('Channel: Es/N0 = %.1f dB | sample SNR = %.2f dB\n', ...
    target_snr_db, input_snr);

%% 5. Receiver
rx = Receiver(M, Rs, fc, Vd, Ns, coefficients.Num);
rx.do_demodulation(r_t, t_axis);
rx.do_lowpass_filter();
rx.do_downsampling(t_axis);
rx.do_matched_filtering();
rx.do_detection();

%% 6. Preamble synchronization and payload error measurement
rx_seq = rx.detection.symbols;
delay_start = finddelay(preamble, rx_seq);
idx_payload_start = delay_start + preamble_len + 1;
idx_payload_end = idx_payload_start + payload_len - 1;

rx.plot_channel_output(r_t, t_axis);
rx.plot_mixer_output();
rx.plot_lpf_and_downsampler();
rx.plot_matched_filter_output();
rx.plot_raw_detected_voltages(preamble_len);

if idx_payload_start > 0 && idx_payload_end <= length(rx_seq)
    rx_payload = rx_seq(idx_payload_start:idx_payload_end);
    errors = sum(payload ~= rx_payload);
    SER = errors / payload_len;
    rx.plot_aligned_payload(idx_payload_start);
    fprintf('Preamble lag: %d symbols\n', delay_start);
    fprintf('Payload SER: %.5f (%d errors / %d symbols)\n', ...
        SER, errors, payload_len);
    fprintf('First 15 transmitted symbols: %s\n', mat2str(payload(1:15)));
    fprintf('First 15 received symbols:    %s\n', mat2str(rx_payload(1:15)));
else
    warning('CommuSystem:SynchronizationFailed', ...
        'The synchronized payload is outside the received frame.');
end

%% 7. Matched-filter eye diagram
safe_start = 100;
safe_end = min(12000, length(rx.matched.signal));
eyediagram(rx.matched.signal(safe_start:safe_end), rx.Ns);
title('Receiver Matched Filter Output');
