classdef Transmitter < handle
    % TRANSMITTER Implements the digital transmitter chain for M-ASK modulation.
    %
    % The Transmitter class encapsulates the entire signal generation process:
    % 1. Random Symbol Generation (Source)
    % 2. Symbol Mapping (Bits -> Voltage Levels)
    % 3. Upsampling (Zero Insertion)
    % 4. Pulse Shaping (RRC / Rectangular Filter)
    % 5. Digital-to-Analog Conversion (DAC Interpolation)
    % 6. RF Modulation (Carrier Mixing)
    %
    % Usage:
    %   tx = Transmitter(M, Rs, fc, Vd, 'RRC');
    %   tx.generate_symbols(1000);
    %   tx.do_mapping();
    %   ...
    %
    % See also: Receiver, run_lab

    %% Properties
    properties
        % System Parameters
        M       % Modulation Order (e.g., 2, 4, 8)
        Rs      % Symbol Rate (Baud)
        Ts      % Symbol Period (1/Rs)
        fc      % Carrier Frequency (Hz)
        Vd      % Voltage Distance between symbols (Volts)
        
        % Data Storage
        symbols % Integer symbol sequence (0 to M-1)
        Ns      % Samples per symbol (Oversampling factor)
        
        % Pulse Shaping Configuration
        pulse_type % Filter type: 'RRC' or 'Rectangular'
        rolloff = 0.22; % Roll-off factor for RRC filter (Beta)
        span = 40;      % Filter span in symbols (Truncation length)
        
        % Signal Stages (Structs containing .signal and .time)
        mapping         % Output of the Mapper
        upsampling      % Output of the Upsampler
        pulse_shaping   % Output of the Pulse Shaping Filter
        dac             % Output of the DAC (High sampling rate)
        modulation      % Final RF Output
    end
    
    %% Public Methods
    methods
        function obj = Transmitter(M, Rs, fc, Vd, pulse_type)
            % TRANSMITTER Constructor
            % Initializes the transmitter with system parameters.
            %
            % Args:
            %   M (int): Modulation order.
            %   Rs (double): Symbol rate in Hz.
            %   fc (double): Carrier frequency in Hz.
            %   Vd (double): Voltage distance between constellation points.
            %   pulse_type (string): 'RRC' (default) or 'Rectangular'.
            
            obj.M = M; 
            obj.Rs = Rs; 
            obj.Ts = 1/Rs; 
            obj.fc = fc; 
            obj.Vd = Vd;
            
            if nargin < 5
                obj.pulse_type = 'RRC'; 
            else
                obj.pulse_type = pulse_type; 
            end
            
            % Initialize Output Structures
            obj.mapping = struct('signal', [], 'time', [], 'constellation', []);
            obj.upsampling = struct('signal', [], 'time', []);
            obj.pulse_shaping = struct('signal', [], 'time', []);
            obj.dac = struct('signal', [], 'time', []);
            obj.modulation = struct('signal', [], 'time', []);
        end
        
        function generate_symbols(obj, NbSymb)
            % GENERATE_SYMBOLS Generates a random sequence of M-ary symbols.
            %
            % Args:
            %   NbSymb (int): Total number of symbols to generate.
            obj.symbols = randi([0, obj.M-1], 1, NbSymb);
        end
        
        function do_mapping(obj)
            % DO_MAPPING Maps integer symbols to voltage levels (PAM/ASK).
            % The mapping is centered around zero.
            % Formula: Level = (2*i - (M-1)) * (Vd/2)
            
            sym_indices = 0:(obj.M-1);
            levels = (2*sym_indices - (obj.M-1)) * (obj.Vd / 2);
            
            obj.mapping.constellation = levels;
            obj.mapping.signal = levels(obj.symbols + 1);
            obj.mapping.time = (0:length(obj.symbols)-1) * obj.Ts;
        end
        
        function do_upsampling(obj, Ns)
            % DO_UPSAMPLING Inserts Ns-1 zeros between each symbol.
            %
            % Args:
            %   Ns (int): Oversampling factor (Samples per Symbol).
            
            obj.Ns = Ns; 
            obj.upsampling.signal = upsample(obj.mapping.signal, Ns);
            
            dt = obj.Ts / Ns;
            obj.upsampling.time = (0:length(obj.upsampling.signal)-1) * dt;
        end
        
        function do_pulse_shaping(obj)
            % DO_PULSE_SHAPING Convolves the upsampled signal with a shaping filter.
            % Supports 'Rectangular' (NRZ) and 'RRC' (Root Raised Cosine).
            
            % 1. Pulse Selection
            if strcmpi(obj.pulse_type, 'Rectangular')
                % Rectangular pulse of duration Ts
                pulse = ones(1, obj.Ns); 
            else
                % Root Raised Cosine Filter
                pulse = rcosdesign(obj.rolloff, obj.span, obj.Ns, 'sqrt'); 
            end
            
            % 2. Convolution
            obj.pulse_shaping.signal = conv(obj.upsampling.signal, pulse);
            
            % 3. Time Axis Update
            dt = obj.Ts / obj.Ns;
            obj.pulse_shaping.time = (0:length(obj.pulse_shaping.signal)-1) * dt;
        end
        
        function do_dac(obj, Fs_target)
            % DO_DAC Simulates Digital-to-Analog conversion via interpolation.
            %
            % Args:
            %   Fs_target (double): Target analog sampling frequency (Hz).
            %                       Must satisfy Nyquist: Fs_target > 2*fc.
            
            t_old = obj.pulse_shaping.time;
            y_old = obj.pulse_shaping.signal;
            
            dt_new = 1/Fs_target;
            t_new = 0:dt_new:max(t_old);
            
            % Use Linear interpolation for simple approximation of an analog signal
            obj.dac.signal = interp1(t_old, y_old, t_new, 'linear', 'extrap');
            obj.dac.time = t_new;
        end
        
        function do_modulation(obj)
            % DO_MODULATION Mixes the baseband signal with a cosine carrier.
            % y(t) = x(t) * cos(2*pi*fc*t)
            
            carrier = cos(2 * pi * obj.fc * obj.dac.time);
            obj.modulation.signal = obj.dac.signal .* carrier;
            obj.modulation.time = obj.dac.time;
        end
        
        %% Visualization Methods
        
        function plot_constellation(obj)
            % PLOT_CONSTELLATION Plots the symbol mapping constellation.
            % Useful for verifying the voltage levels.
            
            if isempty(obj.mapping.constellation)
                warning('Constellation data is empty. Run do_mapping first.');
                return;
            end
            
            h = scatterplot(obj.mapping.constellation);
            grid on;
            
            % Adjust axis limits for better visibility
            max_val = max(abs(obj.mapping.constellation));
            limit = max_val + 1;
            xlim([-limit, limit]); ylim([-limit, limit]);
            
            title(['Tx Constellation (M=' num2str(obj.M) ')']);
        end
        
        function plot_baseband_signals(obj, num_symbols)
            % PLOT_BASEBAND_SIGNALS Superimposes Mapper, Upsampler, and Shaper outputs.
            % Visualizes the transition from discrete symbols to continuous waveforms.
            %
            % Args:
            %   num_symbols (int, optional): Number of symbols to display (default: 20).
            
            if nargin < 2, num_symbols = 20; end 
            
            % Bounds check
            num_symbols = min(num_symbols, length(obj.symbols));
            num_samples = num_symbols * obj.Ns;
            
            figure('Name', 'Tx: Baseband Superimposed'); hold on;
            box on; grid on;
            
            % 1. Mapper Output (Blue circles + stems)
            stem(obj.mapping.time(1:num_symbols), obj.mapping.signal(1:num_symbols), ...
                'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 8, ...
                'LineWidth', 1.5, 'Color', '#0072BD', ... 
                'DisplayName', 'Mapper Output (Symbols)');
            
            % 2. Upsampler Output (Red dashed stems)
            stem(obj.upsampling.time(1:num_samples), obj.upsampling.signal(1:num_samples), ...
                'LineStyle', '--', 'Marker', 'none', ...
                'LineWidth', 1.5, 'Color', '#D95319', ... 
                'DisplayName', 'Upsampler Output (Zeros Inserted)');
            
            % 3. Pulse Shaping Output (Yellow solid line)
            plot(obj.pulse_shaping.time(1:num_samples), obj.pulse_shaping.signal(1:num_samples), ...
                'LineStyle', '-', 'LineWidth', 2, 'Color', '#EDB120', ... 
                'DisplayName', ['Pulse Shaping Filter (' obj.pulse_type ')']);
            
            hold off;
            xlabel('Time (s)'); ylabel('Voltage (V)');
            title('Baseband Processing Stages (Superimposed)');
            legend('Location', 'best');
        end
        
        function plot_modulator_output(obj, time_limit)
            % PLOT_MODULATOR_OUTPUT Compares the Baseband Envelope with the RF Carrier.
            %
            % Args:
            %   time_limit (double, optional): Time duration to plot in seconds.
            
            if nargin < 2, time_limit = 20 * obj.Ts; end 
            
            % --- Data Slicing (Performance Optimization) ---
            idx_end = find(obj.modulation.time >= time_limit, 1);
            if isempty(idx_end), idx_end = length(obj.modulation.time); end
            
            t_slice = obj.modulation.time(1:idx_end);
            dac_slice = obj.dac.signal(1:idx_end);
            mod_slice = obj.modulation.signal(1:idx_end);
            
            figure('Name', 'Tx: Modulator Output Superimposed'); hold on;
            box on; grid on;
            
            % 1. Modulator Output (RF Carrier - Background)
            plot(t_slice, mod_slice, ...
                'Color', [1, 0.7, 0.7], ... % Light Red/Pink
                'LineWidth', 0.5, ...       
                'DisplayName', ['Modulator Output (RF @ ' num2str(obj.fc/1e6) ' MHz)']);
            
            % 2. DAC Output (Baseband Envelope - Foreground)
            plot(t_slice, dac_slice, ...
                'Color', '#0072BD', ...     % Dark Blue
                'LineWidth', 2, ...         
                'DisplayName', 'DAC Output (Baseband Envelope)');
            
            hold off;
            xlabel('Time (s)'); ylabel('Voltage (V)');
            title('Modulation Comparison: Baseband Envelope vs RF Carrier');
            legend('Location', 'best');
            xlim([0, time_limit]); 
        end
    end
end