classdef Receiver < handle
    % RECEIVER Implements the digital receiver chain for M-ASK modulation.
    %
    % The Receiver class encapsulates the signal processing stages to recover
    % symbols from the received RF signal:
    % 1. Coherent Demodulation (Mixing with Local Oscillator)
    % 2. Lowpass Filtering (Removal of 2*fc component and noise)
    % 3. Downsampling (Resampling to baseband rate Rs*Ns)
    % 4. Matched Filtering (RRC Filter to maximize SNR)
    % 5. Symbol Detection (Sampling, AGC, and Voltage Decision)
    %
    % Usage:
    %   rx = Receiver(M, Rs, fc, Vd, Ns, b, a);
    %   rx.do_demodulation(r_noisy, t_axis);
    %   rx.do_lowpass_filter();
    %   ...
    %
    % See also: Transmitter, run_lab

    %% Properties
    properties
        % System Parameters
        M       % Modulation Order
        Rs      % Symbol Rate (Baud)
        Ts      % Symbol Period (1/Rs)
        fc      % Carrier Frequency (Hz)
        Vd      % Voltage Distance (Volts)
        Ns      % Samples per symbol (Baseband oversampling factor)
        
        % Filter Configuration
        rolloff = 0.22; % RRC Roll-off factor
        span = 40;      % RRC Filter span (symbols)
        
        b, a            % Lowpass Filter Coefficients (IIR or FIR)
        mf_coeffs       % Matched Filter Coefficients (RRC)
        
        % Signal Stages (Structs containing .signal and .time)
        demodulation    % Output of the Demodulator (Mixer)
        filtered        % Output of the Lowpass Filter
        downsampled     % Output of the Resampler
        matched         % Output of the Matched Filter
        detection       % Final Detected Symbols and Voltages
    end
    
    %% Public Methods
    methods
        function obj = Receiver(M, Rs, fc, Vd, Ns, b, a)
            % RECEIVER Constructor
            % Initializes the receiver parameters and filters.
            %
            % Args:
            %   M (int): Modulation order.
            %   Rs (double): Symbol rate.
            %   fc (double): Carrier frequency.
            %   Vd (double): Voltage distance.
            %   Ns (int): Samples per symbol for baseband processing.
            %   b, a (arrays): Filter coefficients (numerator, denominator).
            
            obj.M = M; 
            obj.Rs = Rs; 
            obj.fc = fc; 
            obj.Vd = Vd; 
            obj.Ns = Ns;
            obj.Ts = 1/Rs;
            
            % Filter Configuration
            obj.b = b;
            if nargin < 7 || isempty(a)
                obj.a = 1; 
                fprintf('Receiver: Configured for FIR Filter (a=1).\n');
            else
                obj.a = a;
                fprintf('Receiver: Configured for IIR Filter (b, a).\n');
            end
            
            % Generate Matched Filter (RRC)
            obj.mf_coeffs = rcosdesign(obj.rolloff, obj.span, obj.Ns, 'sqrt');
            
            % Initialize Output Structures
            obj.demodulation = struct('signal', [], 'time', []);
            obj.filtered     = struct('signal', [], 'time', []);
            obj.downsampled  = struct('signal', [], 'time', []);
            obj.matched      = struct('signal', [], 'time', []);
            obj.detection    = struct('voltages', [], 'symbols', [], 'phase_idx', []);
        end
        
        function do_demodulation(obj, r_t, t_axis)
            % DO_DEMODULATION Multiplies received signal with local oscillator.
            % y(t) = r(t) * cos(2*pi*fc*t)
            %
            % Args:
            %   r_t (array): Received RF signal (with noise).
            %   t_axis (array): Time vector corresponding to r_t.
            
            lo = cos(2 * pi * obj.fc * t_axis);
            
            % Store in struct
            obj.demodulation.signal = r_t .* lo;
            obj.demodulation.time = t_axis;
        end
        
        function do_lowpass_filter(obj)
            % DO_LOWPASS_FILTER Removes high-frequency components (2*fc).
            % Also multiplies by 2 to recover original baseband amplitude.
            
            if isempty(obj.demodulation.signal)
                error('Receiver: No demodulated signal. Run do_demodulation first.');
            end
            
            % Apply Filter and Gain Correction
            % Note: Coherent demodulation introduces a factor of 0.5, so we multiply by 2.
            obj.filtered.signal = filter(obj.b, obj.a, obj.demodulation.signal) * 2;
            obj.filtered.time = obj.demodulation.time;
        end
        
        function do_downsampling(obj, t_axis_in)
            % DO_DOWNSAMPLING Resamples the signal to baseband rate (Rs * Ns).
            % This handles the case where Fs_analog is not an integer multiple of Rs.
            %
            % Args:
            %   t_axis_in (array): The time vector of the input signal (usually high Fs).
            
            % 1. Calculate Target Sampling Parameters
            Fs_target = obj.Rs * obj.Ns; 
            dt_target = 1 / Fs_target;
            
            % 2. Generate New Time Grid
            t_new = t_axis_in(1) : dt_target : t_axis_in(end);
            
            % 3. Perform Resampling (Interpolation)
            % Accessing .signal from the struct
            sig_resampled = interp1(t_axis_in, obj.filtered.signal, t_new, 'spline');
            
            % 4. Store Result in Struct
            obj.downsampled.signal = sig_resampled;
            obj.downsampled.time = t_new;
        end        
        
        function do_matched_filtering(obj)
            % DO_MATCHED_FILTERING Convolves with the RRC Matched Filter.
            % Maximizes SNR at the optimal sampling instant.
            
            if isempty(obj.downsampled.signal)
                error('Receiver: No downsampled signal. Run do_downsampling first.');
            end
            
            % Convolution
            obj.matched.signal = conv(obj.downsampled.signal, obj.mf_coeffs);
            
            % Update Time Axis (Convolution increases length)
            dt = obj.Ts / obj.Ns;
            obj.matched.time = (0:length(obj.matched.signal)-1) * dt;
        end
        
        function do_detection(obj)
            % DO_DETECTION Performs Blind Phase Sync, AGC, and Symbol Decision.
            % 1. Finds the optimal sampling phase (Max Energy).
            % 2. Scales amplitude (AGC) to match ideal levels.
            % 3. Maps voltages to nearest valid symbols.
            
            sig = obj.matched.signal;
            
            % --- 1. Blind Phase Synchronization (Energy Maximization) ---
            best_phase = 1;
            max_energy = -1;
            
            for p = 1:obj.Ns
                % Subsample the signal at phase p
                temp = sig(p : obj.Ns : end);
                current_energy = sum(abs(temp).^2);
                if current_energy > max_energy
                    max_energy = current_energy;
                    best_phase = p;
                end
            end
            
            % Extract samples at optimal phase
            raw_samples = sig(best_phase : obj.Ns : end);
            
            % --- 2. Blind AGC (Automatic Gain Control) ---
            rx_power = mean(raw_samples.^2);
            
            sym_indices = 0:(obj.M-1);
            ideal_levels = (2*sym_indices - (obj.M-1)) * (obj.Vd/2);
            ideal_power = mean(ideal_levels.^2);
            
            if rx_power > 0
                scale_factor = sqrt(ideal_power / rx_power);
            else
                scale_factor = 1;
            end
            
            samples = raw_samples * scale_factor;
            
            % --- 3. Hard Decision (Nearest Neighbor) ---
            est_syms = zeros(1, length(samples));
            est_volts = zeros(1, length(samples));
            
            for i = 1:length(samples)
                [~, idx] = min(abs(samples(i) - ideal_levels));
                est_syms(i) = sym_indices(idx);
                est_volts(i) = ideal_levels(idx);
            end
            
            % Store results in Struct
            obj.detection.voltages = est_volts;
            obj.detection.symbols = est_syms;
            obj.detection.phase_idx = best_phase;
        end
        
        %% Visualization Methods
        function plot_channel_output(obj, r_rf, t_axis_rf, num_samples)
            % PLOT_CHANNEL_OUTPUT Plots the raw signal received from the AWGN channel.
            % This visualizes the high-frequency RF carrier with noise.
            %
            % Usage:
            %   rx.plot_channel_output(r_t, t_axis, 1000);
            %
            % Args:
            %   r_rf (array): The raw received RF signal (output of awgn).
            %   t_axis_rf (array): The time vector corresponding to r_rf.
            %   num_samples (int, optional): Number of samples to plot (default: 1000).
            
            if nargin < 4 || isempty(num_samples), num_samples = 1000; end
            
            % Slice data from the middle to avoid transients
            [t_slice, y_slice] = obj.get_slice(t_axis_rf, r_rf, num_samples);
            
            figure('Name', 'Rx: 1. AWGN Channel Output (RF)');
            plot(t_slice, y_slice, 'Color', [0.5 0.5 0.5]); % Grey color
            title('AWGN Channel Output (RF Input)');
            ylabel('Amplitude (V)'); xlabel('Time (s)');
            legend(['RF Signal @ ' num2str(obj.fc/1e9) ' GHz'], 'Location', 'best');
            grid on;
        end

        function plot_mixer_output(obj, num_samples)
            % PLOT_MIXER_OUTPUT Plots the signal after the mixing stage.
            % Visualizes the combination of Baseband and 2*fc components.
            %
            % Usage:
            %   rx.plot_mixer_output(1000);
            
            if isempty(obj.demodulation.signal)
                warning('Receiver: No demodulation data found.'); return;
            end
            if nargin < 2 || isempty(num_samples), num_samples = 1000; end
            
            % Slice data
            [t_slice, y_slice] = obj.get_slice(obj.demodulation.time, ...
                                               obj.demodulation.signal, num_samples);
            
            figure('Name', 'Rx: 2. Mixer Output');
            plot(t_slice, y_slice, 'b');
            title('Passband Demodulator Output (Mixer)');
            ylabel('Amplitude (V)'); xlabel('Time (s)');
            legend('Baseband + 2f_c Component', 'Location', 'best');
            grid on;
        end

        function plot_lpf_and_downsampler(obj, num_samples)
            % PLOT_LPF_AND_DOWNSAMPLER Plots the Clean Baseband signal.
            % Visualizes the Low Pass Filter output overlaid with Downsampled points.
            %
            % Usage:
            %   rx.plot_lpf_and_downsampler(1000);
            
            if isempty(obj.filtered.signal)
                warning('Receiver: No filtered data found.'); return;
            end
            if nargin < 2 || isempty(num_samples), num_samples = 1000; end
            
            % 1. Slice the LPF signal (Continuous-like)
            [t_lpf, y_lpf] = obj.get_slice(obj.filtered.time, ...
                                           obj.filtered.signal, num_samples);
            
            % 2. Find corresponding Downsampled points within this time window
            t_start = t_lpf(1);
            t_end = t_lpf(end);
            idx_ds = find(obj.downsampled.time >= t_start & obj.downsampled.time <= t_end);
            
            figure('Name', 'Rx: 3. LPF & Downsampler');
            hold on;
            % Plot LPF Line
            plot(t_lpf, y_lpf, 'k-', 'LineWidth', 1.5, ...
                'DisplayName', 'LPF Output (Clean Baseband)');
            
            % Plot Downsampled Stems
            if ~isempty(idx_ds)
                stem(obj.downsampled.time(idx_ds), obj.downsampled.signal(idx_ds), ...
                    'r', 'Filled', 'LineWidth', 1.2, 'MarkerSize', 6, ...
                    'DisplayName', ['Downsampled (@ ' num2str(obj.Rs/1e6) ' MBaud)']);
            end
            hold off;
            
            title('Low Pass Filter & Downsampler');
            ylabel('Amplitude (V)'); xlabel('Time (s)');
            legend('Location', 'best');
            grid on;
        end

        function plot_matched_filter_output(obj, num_samples)
            % PLOT_MATCHED_FILTER_OUTPUT Plots the waveform after matched filtering.
            % Now uses Time axis for consistency with other stages.
            %
            % Usage:
            %   rx.plot_matched_filter_output(5000);
            
            if nargin < 2 || isempty(num_samples)
                num_samples = 5000;
            end
            
            % Check if data exists
            if isempty(obj.matched.signal)
                warning('Receiver: No matched filter data found.'); return;
            end

            [t_slice, y_slice] = obj.get_slice(obj.matched.time, ...
                                               obj.matched.signal, num_samples);
            
            figure('Name', 'Rx: 4. Matched Filter Waveform'); 
            plot(t_slice, y_slice, 'b-'); 
            
            title('Matched Filter Output (Time Domain)'); 
            xlabel('Time (s)'); 
            ylabel('Amplitude');
            grid on;
            
            % Optional: Add a legend to show sampling info
            legend(['Signal @ ' num2str(obj.Rs * obj.Ns / 1e6) ' MHz'], 'Location', 'best');
        end
        
        function plot_raw_detected_voltages(obj, num_symbols)
            % PLOT_RAW_DETECTED_VOLTAGES Plots the detected voltages (Preamble + Payload).
            % Args:
            %   num_symbols (int, optional): Number of symbols to plot.
            
            if nargin < 2 || isempty(num_symbols)
                num_symbols = 50;
            end
            
            limit = min(num_symbols, length(obj.detection.voltages));
            
            figure('Name', 'Rx: Raw Detected Voltages'); 
            stem(obj.detection.voltages(1:obj.span + 1 + limit)); 
            xlim([obj.span + 1, obj.span + 1 + limit]);
            title('Raw Detected Voltages (Preamble + Payload)'); 
            xlabel('Symbol Index (Raw)'); 
            ylabel('Voltage Level');
            grid on;
        end
        
        function plot_aligned_payload(obj, payload_start_idx, num_symbols)
            % PLOT_ALIGNED_PAYLOAD Plots only the payload part of the signal.
            % Args:
            %   payload_start_idx (int): Index where payload starts.
            %   num_symbols (int, optional): Number of symbols to plot.
            
            if nargin < 2 || isempty(payload_start_idx)
                error('Receiver: plot_aligned_payload requires payload_start_idx.');
            end
            
            if nargin < 3 || isempty(num_symbols)
                num_symbols = 50;
            end
            
            % Safety check
            if payload_start_idx < length(obj.detection.voltages)
                
                % Accessing .voltages from struct
                payload_volts = obj.detection.voltages(payload_start_idx : end);
                limit = min(num_symbols, length(payload_volts));
                
                figure('Name', 'Rx: Payload Only (Aligned)');
                stem(payload_volts(1 : limit), 'Filled', 'LineWidth', 1.5);
                
                title('Detected Payload Voltages (Preamble Removed)');
                ylabel('Voltage Level');
                xlabel('Payload Symbol Index');
                grid on;
            else
                warning('Receiver: payload_start_idx is out of bounds. Skipping Payload plot.');
            end
        end
        
        %% Helper Method (Private)
        function [t_slice, y_slice] = get_slice(~, t_full, y_full, n_points)
            % GET_SLICE Helper to extract a middle segment of data.
            total = length(t_full);
            start_idx = floor(total/2);
            end_idx = min(start_idx + n_points, total);
            range = start_idx : end_idx;
            
            t_slice = t_full(range);
            y_slice = y_full(range);
        end
    end
end