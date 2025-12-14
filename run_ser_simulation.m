%% Monte Carlo Simulation for SER vs SNR (Preamble Sync Version)
% This script simulates the SER using a Robust Barker Preamble for synchronization,
% matching the realistic logic used in the lab experiments.

clear; close all; clc;
load("my_filter_coeffs_FIR.mat"); % 确保工作区有 Num 变量

%% 1. System Configuration
M = 4;                  % Modulation Order (4-ASK)
Rs = 3.84e6;               % Symbol Rate
fc = 1920e6;               % Carrier Frequency
Vd = 1;                 % Voltage Distance
Ns = 64;                % Samples per symbol

% 采样率设置
% Fs_baseband = Rs * Ns;    
Fs_analog = 4500e6; 

% Simulation Parameters
SNR_dB_range = 0:1:20;    % SNR range
payload_len = 100000;     % 只统计 Payload 的误码率
                          % 注意：实际发送的符号数 = preamble_len + payload_len

% Initialize Results
SER_results = zeros(size(SNR_dB_range));

%% 2. Simulation Loop
hWait = waitbar(0, 'Running Monte Carlo Simulation...');

for i = 1:length(SNR_dB_range)
    current_snr = SNR_dB_range(i);
    
    % Update Waitbar
    waitbar(i/length(SNR_dB_range), hWait, ...
        sprintf('Simulating SNR = %d dB...', current_snr));
    
    % --- A. Transmitter ---
    tx = Transmitter(M, Rs, fc, Vd, 'RRC'); 
    
    % [核心修改 1]：构造带有巴克码前导码的帧
    % ------------------------------------------------------
    % 1. 生成巴克码 (Barker-13 x 4)
    barker13 = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
    seq_binary = [barker13, barker13, barker13, barker13]; 
    
    % 2. 映射到 M-ASK 的最大能量点 (0 和 M-1)
    preamble = zeros(1, length(seq_binary));
    preamble(seq_binary == 1)  = M - 1; % 最高电平
    preamble(seq_binary == -1) = 0;     % 最低电平
    
    % 3. 生成随机 Payload
    payload = randi([0 M-1], 1, payload_len);
    
    % 4. 组合成发送帧，并手动注入到 tx 对象
    tx.symbols = [preamble, payload]; 
    % ------------------------------------------------------
    
    % 继续发送流程
    tx.do_mapping();
    tx.do_upsampling(Ns);
    tx.do_pulse_shaping();
    tx.do_dac(Fs_analog);
    tx.do_modulation();
    
    % --- B. AWGN Channel ---
    % 计算处理增益并修正 SNR
    OversamplingFactor = Fs_analog / tx.Rs;
    process_gain = 10 * log10(OversamplingFactor);
    snr_corrected = current_snr - process_gain;
    
    r_noisy = awgn(tx.modulation.signal, snr_corrected, 'measured');
    t_axis = tx.modulation.time;
    
    % --- C. Receiver ---
    % 使用加载的 FIR 系数 Num
    rx = Receiver(M, Rs, fc, Vd, Ns, Num);
    
    % Process
    rx.do_demodulation(r_noisy, t_axis);
    rx.do_lowpass_filter();
    rx.do_downsampling(t_axis);
    rx.do_matched_filtering();
    rx.do_detection();
    
    % --- D. Synchronization (Using Preamble) ---
    % [核心修改 2]：只使用前导码进行同步，模拟真实接收机逻辑
    
    % 1. 在接收到的符号流中寻找前导码的位置
    % finddelay 内部本质是互相关，它会告诉我们 rx_syms 滞后了 preamble 多少个点
    rx_detected_all = rx.detected_symbols;
    
    % delay_idx 是前导码在接收序列中的起始索引 (0-based 还是 1-based 取决于实现，MATLAB finddelay 返回延迟)
    % 如果 rx 从第 10 个点开始是 preamble，finddelay 返回 9 (因为 delay 9 个点)
    % 所以起始位置通常是 delay + 1
    delay = finddelay(preamble, rx_detected_all);
    
    % 2. 提取 Payload
    pre_len = length(preamble);
    pay_len = length(payload);
    
    % 计算 Payload 在接收序列中的预计起始和结束位置
    % 逻辑：接收序列起始位置 = delay + 1
    % Payload 起始位置 = 接收序列起始位置 + Preamble长度
    idx_payload_start = delay + 1 + pre_len;
    idx_payload_end   = idx_payload_start + pay_len - 1;
    
    % 3. 安全性检查与误码统计
    % 必须确保索引没有越界 (在低信噪比下同步可能算错导致越界)
    if idx_payload_start > 0 && idx_payload_end <= length(rx_detected_all)
        
        % 提取接收到的 Payload
        rx_payload = rx_detected_all(idx_payload_start : idx_payload_end);
        
        % 计算误码 (只比较 Payload，不比较 Preamble)
        num_errors = sum(payload ~= rx_payload);
        
    else
        % 如果同步完全失败 (找不到头或者索引越界)
        % 这种情况通常发生在 SNR 极低时。我们判为全部错误或丢失帧。
        num_errors = pay_len; % 惩罚：认为所有位都错了
        % 或者 num_errors = round(0.5 * pay_len); % 盲猜也有 50% 错误率
    end
    
    SER = num_errors / pay_len;
    SER_results(i) = SER;
    
    fprintf('SNR: %2d dB | Errors: %4d / %d | SER: %.2e\n', ...
            current_snr, num_errors, pay_len, SER);
end
close(hWait);

%% 3. Plotting (保持不变)
figure('Name', 'SER vs SNR Performance (Preamble Sync)');
SNR_lin = 10.^(SNR_dB_range / 10); 
my_qfunc = @(x) 0.5 * erfc(x / sqrt(2));
theo_SER = 2 * ((M-1)/M) * my_qfunc(sqrt(3 * SNR_lin / (M^2 - 1)));

semilogy(SNR_dB_range, SER_results, 'b-o', 'LineWidth', 2, ...
    'MarkerFaceColor', 'b', 'DisplayName', 'Simulation (Preamble Sync)');
hold on;
semilogy(SNR_dB_range, theo_SER, 'k-', 'LineWidth', 1.5, ...
    'DisplayName', ['Theory (' num2str(M) '-ASK)']);

grid on;
title(['SER for ' num2str(M) '-ASK with Barker Sync']);
xlabel('SNR (E_s/N_0) [dB]');
ylabel('Symbol Error Rate (SER)');
legend('Location', 'southwest');
axis([min(SNR_dB_range) max(SNR_dB_range) 1e-6 1]);
hold off;