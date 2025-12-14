classdef Receiver < handle
    % RECEIVER Class (Final Clean Version)
    % Implements: Demod -> IIR LPF -> Downsample -> Matched Filter -> Blind Detection
    
    properties
        M, Rs, Ts, fc, Vd, Ns
        rolloff = 0.22; span = 40;
        
        b, a        % IIR Filter Coefficients
        mf_coeffs   % Matched Filter (RRC)
        
        % Signals
        demodulated, filtered, downsampled, matched
        detected_symbols, detected_voltages
    end
    
    methods
        function obj = Receiver(M, Rs, fc, Vd, Ns, b, a)
            obj.M = M; obj.Rs = Rs; obj.fc = fc; obj.Vd = Vd; obj.Ns = Ns;
            obj.Ts = 1/Rs;
            obj.b = b;
            if nargin < 7 || isempty(a)
                obj.a = 1; 
                fprintf('Receiver: Configured for FIR Filter (a=1).\n');
            else
                obj.a = a;
                fprintf('Receiver: Configured for IIR Filter (b, a).\n');
            end
            obj.mf_coeffs = rcosdesign(obj.rolloff, obj.span, obj.Ns, 'sqrt');
        end
        
        function do_demodulation(obj, r_t, t_axis)
            lo = cos(2 * pi * obj.fc * t_axis);
            obj.demodulated = r_t .* lo;
        end
        
        function do_lowpass_filter(obj)
            % Apply IIR Filter
            obj.filtered = filter(obj.b, obj.a, obj.demodulated) * 2;
        end
        
        function do_downsampling(obj, t_axis_in)
    % DO_DOWNSAMPLING (Array Version)
    % Performs resampling using interpolation to handle non-integer factors.
    % Adapted for the "Clean Version" where signals are Arrays, not Structs.
    
    % 1. Calculate Target Sampling Parameters
    Fs_target = obj.Rs * obj.Ns; 
    dt_target = 1 / Fs_target;
    
    % 2. Generate New Time Grid
    % Create a time vector for the desired output sample rate
    t_new = t_axis_in(1) : dt_target : t_axis_in(end);
    
    % 3. Perform Resampling (Interpolation)
    % FIX: Use 'obj.filtered' directly (it is an Array, remove .signal)
    sig_resampled = interp1(t_axis_in, obj.filtered, t_new, 'spline');
    
    % 4. Store Result
    % FIX: Store directly into 'obj.downsampled' (as an Array)
    obj.downsampled = sig_resampled;
    
    % Note: In the Array version, we don't store 't_new' in the object.
    % The matched filter and detection steps work on array indices, 
    % so strict time-axis tracking inside the object is not required 
    % for the BER calculation.
end        
        function do_matched_filtering(obj)
            obj.matched = conv(obj.downsampled, obj.mf_coeffs);
        end
        
        function do_detection(obj)
            sig = obj.matched;
            Ns = obj.Ns;
            
            % --- 1. Blind Phase Synchronization ---
            best_phase = 1;
            max_energy = -1;
            
            for p = 1:Ns
                temp = sig(p : Ns : end);
                current_energy = sum(abs(temp).^2);
                if current_energy > max_energy
                    max_energy = current_energy;
                    best_phase = p;
                end
            end
            
            % Extract samples at optimal phase
            raw_samples = sig(best_phase : Ns : end);
            
            % --- 2. Blind AGC ---
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
            % obj.matched = obj.matched * scale_factor;
            % --- 3. Decision ---
            est_syms = zeros(1, length(samples));
            est_volts = zeros(1, length(samples));
            
            for i = 1:length(samples)
                [~, idx] = min(abs(samples(i) - ideal_levels));
                est_syms(i) = sym_indices(idx);
                est_volts(i) = ideal_levels(idx);
            end
            
            obj.detected_symbols = est_syms;
            obj.detected_voltages = est_volts;
        end
        
        function plot_all_separate(obj, payload_start_idx)
    % PLOT_ALL_SEPARATE Plots the internal signals of the receiver.
    %
    % Usage:
    %   rx.plot_all_separate()
    %       Plots the raw matched filter output and detected voltages
    %       (contains the Preamble at the beginning).
    %
    %   rx.plot_all_separate(start_idx)
    %       Additionally plots the "Clean Payload", skipping the Preamble.
    %       This helps verify alignment with the console text output.
    %
    % Args:
    %   payload_start_idx (optional): The integer index where the valid 
    %                                 Payload starts (calculated in run_lab).

    % 1. Plot the Matched Filter Output (Raw Waveform)
    % This signal contains the preamble, payload, and filter transients.
    figure('Name', 'Rx: Matched Filter Waveform (Full Frame)'); 
    plot(obj.matched(1:min(500,end))); 
    title('Matched Filter Output (Includes Preamble)'); 
    xlabel('Sample Index'); ylabel('Amplitude');
    grid on;
    
    % 2. Plot the Raw Detected Voltages (Full Sequence)
    % This visualization starts at index 1, so it shows the Preamble first.
    % e.g., -0.5, -0.5... (if Preamble is 1, 1...)
    figure('Name', 'Rx: Raw Detected Voltages'); 
    stem(obj.detected_voltages(1:min(50,end))); 
    title('Raw Detected Voltages (Preamble + Payload)'); 
    xlabel('Symbol Index (Raw)'); ylabel('Voltage Level');
    grid on;
    
    % 3. (Optional) Plot the Aligned Payload
    % If the caller provided the 'payload_start_idx', we crop the Preamble.
    if nargin > 1 && ~isempty(payload_start_idx)
        
        % Safety check: Ensure the index is within the data bounds
        if payload_start_idx < length(obj.detected_voltages)
            
            % Extract the Payload portion (The "Meat" of the frame)
            payload_volts = obj.detected_voltages(payload_start_idx : end);
            
            figure('Name', 'Rx: Payload Only (Aligned)');
            % Plot the first 50 symbols of the actual Payload.
            % The first value here should match the first value in your Tx log.
            stem(payload_volts(1 : min(50, end)), 'Filled', 'LineWidth', 1.5);
            
            title('Detected Payload Voltages (Preamble Removed)');
            ylabel('Voltage Level (-1.5, -0.5, 0.5, 1.5)');
            xlabel('Payload Symbol Index');
            grid on;
            
        else
            warning('Receiver: payload_start_idx is out of bounds. Skipping Payload plot.');
        end
    end
end
    end
end