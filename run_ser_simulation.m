%% Monte Carlo Simulation for SER vs. SNR (Preamble Synchronization Version)
% This script evaluates the Symbol Error Rate (SER) performance of the M-ASK 
% communication system using Monte Carlo simulations.
%
% Key Features:
% 1. Realistic Frame Structure: Uses a specific Preamble + Payload structure.
% 2. Robust Synchronization: Utilizes a Barker sequence mapped to extreme 
%    voltage levels for robust correlation-based detection.
% 3. Physics-based Channel Modeling: Corrects for oversampling processing gain 
%    to ensure accurate Es/N0 simulation.
% 4. Benchmarking: Compares simulation results against the theoretical M-ASK curve.

clear; close all; clc;
rng(42, 'twister');     % Reproducible payloads and channel noise

% Load pre-designed FIR filter coefficients (variable 'Num')
load("my_filter_coeffs_FIR.mat"); 

%% 1. System Configuration
M = 16;                % Modulation order (16-ASK); try 2, 4, or 8 as well
Rs = 3.84e6;            % Symbol Rate (3.84 MBaud)
fc = 1920e6;            % Carrier Frequency (1.92 GHz)
Vd = 1;                 % Voltage Distance between symbols
Ns = 64;                % Samples per symbol (Baseband Oversampling)

% Sampling Rate Configuration  
Fs_analog = 4500e6;     % DAC/RF Sampling Frequency (4.5 GHz)

% Simulation Parameters
SNR_dB_range = 0:1:35;  % SNR range in dB (Es/N0)
payload_len = 2000;     % Symbols per SNR point; increase for finer SER estimates
                        % Note: Total symbols = preamble_len + payload_len

% Initialize Results Storage
SER_results = zeros(size(SNR_dB_range));

%% 2. Simulation Loop
hWait = waitbar(0, 'Running Monte Carlo Simulation...');

for i = 1:length(SNR_dB_range)
    current_snr = SNR_dB_range(i);
    
    % Update Progress Bar
    waitbar(i/length(SNR_dB_range), hWait, ...
        sprintf('Simulating SNR = %d dB...', current_snr));
    
    % --- A. Transmitter Setup ---
    tx = Transmitter(M, Rs, fc, Vd, 'RRC'); 
    
    % --- Frame Construction with Robust Preamble ---
    % ------------------------------------------------------
    % 1. Generate Barker Sequence (Barker-13 repeated 4 times)
    % Repeated sequences improve correlation peak detection.
    barker13 = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
    seq_binary = [barker13, barker13, barker13, barker13]; 
    
    % 2. Map Preamble to Maximum Energy Levels
    % To maximize detection probability, map bits to the outer-most symbols:
    % Bit 1 -> Symbol (M-1) [Highest Voltage]
    % Bit 0 -> Symbol 0     [Lowest Voltage]
    preamble = zeros(1, length(seq_binary));
    preamble(seq_binary == 1)  = M - 1; 
    preamble(seq_binary == -1) = 0;     
    
    % 3. Generate Random Payload
    payload = randi([0 M-1], 1, payload_len);
    
    % 4. Assemble Frame and Inject into Transmitter
    tx.symbols = [preamble, payload]; 
    % ------------------------------------------------------
    
    % Execute Transmitter Processing Chain
    tx.do_mapping();
    tx.do_upsampling(Ns);
    tx.do_pulse_shaping();
    tx.do_dac(Fs_analog);
    tx.do_modulation();
    
    % --- B. AWGN Channel & Noise Calibration ---
    % The 'awgn' function adds noise over the entire sampling bandwidth (Fs_analog).
    % However, the receiver's matched filter limits the noise bandwidth to approx Rs.
    % For a real RF waveform: Es/N0 = (Fs_analog / (2*Rs)) * SNR_sample.
    % The factor of two follows the real-noise PSD convention N0/2.
    
    OversamplingFactor = Fs_analog / tx.Rs;
    process_gain = 10 * log10(OversamplingFactor / 2);
    snr_corrected = current_snr - process_gain;
    
    r_noisy = awgn(tx.modulation.signal, snr_corrected, 'measured');
    t_axis = tx.modulation.time;
    
    % --- C. Receiver Processing ---
    % Initialize Receiver with the loaded FIR coefficients (Num)
    rx = Receiver(M, Rs, fc, Vd, Ns, Num);
    
    % Execute Receiver Processing Chain
    rx.do_demodulation(r_noisy, t_axis);
    rx.do_lowpass_filter();
    rx.do_downsampling(t_axis);
    rx.do_matched_filtering();
    rx.do_detection();
    
    % --- D. Preamble-Based Synchronization ---
    % Simulates a realistic receiver that knows only the Preamble structure.
    
    % 1. Locate Preamble in the received stream
    % 'finddelay' uses cross-correlation to find the lag between the known
    % preamble and the received symbol stream.
    rx_detected_all = rx.detection.symbols; 
    
    delay = finddelay(preamble, rx_detected_all);
    
    % 2. Calculate Payload Boundaries
    % Logic: Start index = Delay lag + 1 + Length of Preamble
    pre_len = length(preamble);
    pay_len = length(payload);
    
    idx_payload_start = delay + 1 + pre_len;
    idx_payload_end   = idx_payload_start + pay_len - 1;
    
    % 3. Boundary Check and Error Counting
    % Ensure synchronization was successful and indices are valid.
    if idx_payload_start > 0 && idx_payload_end <= length(rx_detected_all)
        
        % Extract the detected payload
        rx_payload = rx_detected_all(idx_payload_start : idx_payload_end);
        
        % Count Symbol Errors
        num_errors = sum(payload ~= rx_payload);
        
    else
        % Synchronization Failure Handling
        % If SNR is too low, correlation might fail, returning invalid indices.
        % We penalize this by marking the entire frame as errors.
        num_errors = pay_len; 
    end
    
    % Calculate and Store SER
    SER = num_errors / pay_len;
    SER_results(i) = SER;
    
    fprintf('SNR: %2d dB | Errors: %4d / %d | SER: %.2e\n', ...
            current_snr, num_errors, pay_len, SER);
end
close(hWait);

%% 3. Visualization and Validation
figure('Name', 'SER vs SNR Performance (Preamble Sync)');

% --- Theoretical Curve Calculation ---
SNR_lin = 10.^(SNR_dB_range / 10); 
% M-ASK Theoretical SER Formula:
% Ps = 2 * ((M-1)/M) * Q( sqrt( 6 * EsN0 / (M^2-1) ) )
% This ideal coherent M-PAM result excludes synchronization and AGC errors.
theo_SER = 2 * ((M-1)/M) * qfunc(sqrt(6 * SNR_lin / (M^2 - 1)));

% --- Plotting ---
% 1. Simulation Results (Blue Markers)
semilogy(SNR_dB_range, SER_results, 'b-o', 'LineWidth', 2, ...
    'MarkerFaceColor', 'b', 'DisplayName', 'Simulation (Preamble Sync)');
hold on;

% 2. Theoretical Limit (Black Line)
semilogy(SNR_dB_range, theo_SER, 'k-', 'LineWidth', 1.5, ...
    'DisplayName', ['Theory (' num2str(M) '-ASK)']);

% 3. Formatting
grid on;
title(['SER Performance for ' num2str(M) '-ASK with Barker Synchronization']);
xlabel('SNR (E_s/N_0) [dB]');
ylabel('Symbol Error Rate (SER)');
legend('Location', 'southwest');
axis([min(SNR_dB_range) max(SNR_dB_range) 1e-6 1]);
hold off;
