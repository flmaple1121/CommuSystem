classdef Transmitter < handle
    % TRANSMITTER Class (Final Standard Version)
    % Implements: Source -> Mapper -> Upsampler -> Pulse Shaping -> DAC -> Modulator
    
    properties
        M, Rs, Ts, fc, Vd   % System Params
        symbols, Ns         % Data
        pulse_type          % 'RRC' or 'Rectangular'
        
        % Filter Params
        rolloff = 0.22;
        span = 40;
        
        % Outputs
        mapping, upsampling, pulse_shaping, dac, modulation
    end

    methods
        function obj = Transmitter(M, Rs, fc, Vd, pulse_type)
            obj.M = M; obj.Rs = Rs; obj.Ts = 1/Rs; obj.fc = fc; obj.Vd = Vd;
            if nargin < 5, obj.pulse_type = 'RRC'; else, obj.pulse_type = pulse_type; end
            
            % Initialize Structs
            obj.mapping = struct('signal', [], 'time', [], 'constellation', []);
            obj.upsampling = struct('signal', [], 'time', []);
            obj.pulse_shaping = struct('signal', [], 'time', []);
            obj.dac = struct('signal', [], 'time', []);
            obj.modulation = struct('signal', [], 'time', []);
        end

        function generate_symbols(obj, NbSymb)
            obj.symbols = randi([0, obj.M-1], 1, NbSymb);
        end

        function do_mapping(obj)
            sym_indices = 0:(obj.M-1);
            levels = (2*sym_indices - (obj.M-1)) * (obj.Vd / 2);
            obj.mapping.constellation = levels;
            obj.mapping.signal = levels(obj.symbols + 1);
            obj.mapping.time = (0:length(obj.symbols)-1) * obj.Ts;
        end

        function do_upsampling(obj, Ns)
            obj.Ns = Ns; 
            obj.upsampling.signal = upsample(obj.mapping.signal, Ns);
            dt = obj.Ts / Ns;
            obj.upsampling.time = (0:length(obj.upsampling.signal)-1) * dt;
        end

        function do_pulse_shaping(obj)
            % Pulse Selection
            if strcmpi(obj.pulse_type, 'Rectangular')
                pulse = ones(1, obj.Ns); % Part A default
            else
                pulse = rcosdesign(obj.rolloff, obj.span, obj.Ns, 'sqrt'); % Part B matched
            end
            
            obj.pulse_shaping.signal = conv(obj.upsampling.signal, pulse);
            dt = obj.Ts / obj.Ns;
            obj.pulse_shaping.time = (0:length(obj.pulse_shaping.signal)-1) * dt;
        end

        function do_dac(obj, Fs_target)
            t_old = obj.pulse_shaping.time;
            y_old = obj.pulse_shaping.signal;
            dt_new = 1/Fs_target;
            t_new = 0:dt_new:max(t_old);
            
            % Use Linear interpolation for better RRC shape preservation
            obj.dac.signal = interp1(t_old, y_old, t_new, 'linear', 'extrap');
            obj.dac.time = t_new;
        end

        function do_modulation(obj)
            carrier = cos(2 * pi * obj.fc * obj.dac.time);
            obj.modulation.signal = obj.dac.signal .* carrier;
            obj.modulation.time = obj.dac.time;
        end
        
% --- Required Plotting Methods (Superimposed Style / English Version) ---
        
        function plot_constellation(obj)
            % PLOT_CONSTELLATION Plots the symbol mapping constellation.
            
            if isempty(obj.mapping.constellation)
                warning('Constellation data is empty. Run do_mapping first.');
                return;
            end
            
            % Use scatterplot for standard visualization
            h = scatterplot(obj.mapping.constellation);
            grid on;
            
            % Adjust axis limits for better visibility
            max_val = max(abs(obj.mapping.constellation));
            limit = max_val + 1;
            xlim([-limit, limit]); ylim([-limit, limit]);
            
            title(['Tx Constellation (M=' num2str(obj.M) ')']);
        end

        function plot_baseband_signals(obj, num_symbols)
            % PLOT_BASEBAND_SIGNALS Plots Mapper, Upsampler, and Shaper on one figure.
            % This visualizes the transition from discrete symbols to continuous waveform.
            
            if nargin < 2, num_symbols = 20; end % Default: show first 20 symbols
            
            % Ensure we don't exceed array bounds
            num_symbols = min(num_symbols, length(obj.symbols));
            num_samples = num_symbols * obj.Ns;
            
            figure('Name', 'Tx: Baseband Superimposed'); hold on;
            box on; grid on;
            
            % 1. Mapper Output (Blue circles + stems)
            % Represents the discrete voltage levels
            stem(obj.mapping.time(1:num_symbols), obj.mapping.signal(1:num_symbols), ...
                'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 8, ...
                'LineWidth', 1.5, 'Color', '#0072BD', ... % Blue
                'DisplayName', 'Mapper Output (Symbols)');
            
            % 2. Upsampler Output (Red dashed stems)
            % Represents the signal after zero-insertion
            stem(obj.upsampling.time(1:num_samples), obj.upsampling.signal(1:num_samples), ...
                'LineStyle', '--', 'Marker', 'none', ...
                'LineWidth', 1.5, 'Color', '#D95319', ... % Red/Orange
                'DisplayName', 'Upsampler Output (Zeros Inserted)');
            
            % 3. Pulse Shaping Output (Yellow solid line)
            % Represents the filtered, continuous baseband signal
            plot(obj.pulse_shaping.time(1:num_samples), obj.pulse_shaping.signal(1:num_samples), ...
                'LineStyle', '-', 'LineWidth', 2, 'Color', '#EDB120', ... % Yellow/Gold
                'DisplayName', ['Pulse Shaping Filter (' obj.pulse_type ')']);
            
            hold off;
            xlabel('Time (s)'); ylabel('Voltage (V)');
            title('Baseband Processing Stages (Superimposed)');
            legend('Location', 'best');
        end

        function plot_modulator_output(obj, time_limit)
            % PLOT_MODULATOR_OUTPUT Superimposes Baseband Envelope and RF Carrier.
            % This visualizes how the baseband signal becomes the envelope of the carrier.
            
            if nargin < 2, time_limit = 20 * obj.Ts; end 
            
            % --- Data Slicing (Performance Optimization) ---
            % Find the index corresponding to the time limit to avoid plotting millions of points
            idx_end = find(obj.modulation.time >= time_limit, 1);
            if isempty(idx_end), idx_end = length(obj.modulation.time); end
            
            % Extract the slice
            t_slice = obj.modulation.time(1:idx_end);
            dac_slice = obj.dac.signal(1:idx_end);
            mod_slice = obj.modulation.signal(1:idx_end);
            
            figure('Name', 'Tx: Modulator Output Superimposed'); hold on;
            box on; grid on;
            
            % 1. Modulator Output (RF Carrier - Background)
            % Plot first as a light-colored background. 
            % The dense high-frequency lines will look like a filled block.
            plot(t_slice, mod_slice, ...
                'Color', [1, 0.7, 0.7], ... % Light Red/Pink
                'LineWidth', 0.5, ...       % Thin lines for density
                'DisplayName', ['Modulator Output (RF @ ' num2str(obj.fc/1e6) ' MHz)']);
            
            % 2. DAC Output (Baseband Signal - Foreground)
            % Plot on top to clearly show it acts as the envelope.
            plot(t_slice, dac_slice, ...
                'Color', '#0072BD', ...     % Dark Blue
                'LineWidth', 2, ...         % Thick line for visibility
                'DisplayName', 'DAC Output (Baseband Envelope)');
            
            hold off;
            xlabel('Time (s)'); ylabel('Voltage (V)');
            title('Modulation Comparison: Baseband Envelope vs RF Carrier');
            legend('Location', 'best');
            xlim([0, time_limit]); % Lock the x-axis view
        end
    end
end