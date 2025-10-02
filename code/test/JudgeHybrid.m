clear;
clc;
close all;

CROSSOVER_N_TO_TEST = 416; %之前的测试结果

fprintf('Testing Hybrid FFT with Crossover_N = %d\n', CROSSOVER_N_TO_TEST);

% 定义测试的FFT规模
% 扩大测试范围以包含你出错的N值
min_pow = 9;
max_pow = 20; % 扩展到 2^20
N_values = 2.^(min_pow:max_pow);

fprintf('This may take a few minutes...\n\n');

% --- 2. 循环测试与数据收集 ---
num_N_values = length(N_values);
results_time = zeros(num_N_values, 4); % 1:Hybrid, 2:Stockham, 3:Radix4-Iter, 4:Built-in
results_error = zeros(num_N_values, 3); % 1:Hybrid, 2:Stockham, 3:Radix4-Iter

% 为了代码清晰，将函数句柄和名称放在这里
base_func_small = @fft_radix4_iterative; 
base_func_large = @stockham_fft;
base_func_small_name = 'Radix-4-Iterative';
base_func_large_name = 'Stockham FFT';


for i = 1:num_N_values
    N = N_values(i);
    fprintf('Testing N = %d (2^%d)...\n', N, log2(N));

    % 生成输入数据
    X_row = randn(1, N) + 1i * randn(1, N);
    X_col = X_row.'; % Stockham FFT 可能需要列向量
    Y_ref = fft(X_row);

    % --- 1. 测试 Hybrid FFT ---
    f_hybrid = @() hybrid_fft(X_row, CROSSOVER_N_TO_TEST);
    [results_time(i, 1), ~, results_error(i, 1)] = run_and_time_fft(f_hybrid, Y_ref);
    


    % --- 3. 测试 Radix-4-Iterative FFT (小规模基础算法) ---
    f_radix4 = @() base_func_small(X_row);
    [results_time(i, 3), ~, results_error(i, 3)] = run_and_time_fft(f_radix4, Y_ref);

    % --- 4. 测试 MATLAB Built-in FFT ---
    f_builtin = @() fft(X_row);
    results_time(i, 4) = timeit(f_builtin);
    
    % 打印单行结果，方便实时观察
    fprintf('  [OK] Finished.\n');
end

% --- 3. 可视化性能对比 ---
figure;
h = loglog(N_values, results_time, '-o', 'LineWidth', 1.5);
grid on;
title(sprintf('Hybrid FFT Performance vs. Base Algorithms (Crossover at N=%d)', CROSSOVER_N_TO_TEST));
xlabel('FFT Size (N)');
ylabel('Execution Time (seconds)');
legend('Hybrid FFT', base_func_large_name, base_func_small_name, 'MATLAB built-in fft', 'Location', 'northwest');
set(h(1), 'LineWidth', 2.5, 'Marker', 's', 'MarkerSize', 8); % 突出显示混合算法的曲线
ylim([min(results_time(:))*0.8, max(results_time(:))*1.2]); % 调整Y轴范围以便观察

% --- 4. 输出数据表格 ---
fprintf('\n================== Final Performance & Correctness Summary ==================\n');
fprintf('%-10s | %-25s | %-25s | %-25s | %-12s\n', 'N', 'Hybrid FFT', base_func_large_name, base_func_small_name, 'Built-in FFT');
fprintf('%-10s | %-12s %-12s | %-12s %-12s | %-12s %-12s | %-12s\n', '', 'Time (s)', 'Error', 'Time (s)', 'Error', 'Time (s)', 'Error', 'Time (s)');
fprintf(repmat('-', 1, 120));
fprintf('\n');

for i = 1:num_N_values
    fprintf('%-10d | %-12.6f %-12.2e | %-12.6f %-12.2e | %-12.6f %-12.2e | %-12.6f\n', ...
        N_values(i), ...
        results_time(i, 1), results_error(i, 1), ...
        results_time(i, 2), results_error(i, 2), ...
        results_time(i, 3), results_error(i, 3), ...
        results_time(i, 4));
end
fprintf(repmat('=', 1, 120));
fprintf('\n\n');

%% --- 辅助函数 (关键修改在此) ---
% 将此函数放在脚本文件末尾，或者保存为 "run_and_time_fft.m"
function [exec_time, result, error_norm] = run_and_time_fft(func_handle, Y_ref, is_col_output)
    % 该函数封装了计时、运行、错误计算和异常处理的逻辑
    %
    % func_handle: 指向要测试的FFT函数的句柄 e.g., @() my_fft(data)
    % Y_ref: 用于计算误差的黄金标准结果 (行向量)
    % is_col_output: (可选) 布尔值，如果函数输出为列向量则设为 true

    if nargin < 3
        is_col_output = false;
    end
    
    norm_Y_ref = norm(Y_ref);

    try
        % ======================= 核心修正 =======================
        % 步骤 1: 先运行一次函数以获取其输出结果，用于正确性检查
        result = func_handle();
        
        % 步骤 2: 使用 timeit 来进行可靠的性能计时
        exec_time = timeit(func_handle);
        % ========================================================

        % 如果函数输出是列向量，转置它以进行误差比较
        if is_col_output
            result = result.';
        end

        % 检查尺寸是否匹配
        if ~isequal(size(result), size(Y_ref))
            warning('Output dimension mismatch for N=%d. Expected [%d %d], got [%d %d].', ...
                    length(Y_ref), size(Y_ref,1), size(Y_ref,2), size(result,1), size(result,2));
            error_norm = Inf; % 标记为错误
        else
            % 计算相对误差
            error_norm = norm(result - Y_ref) / norm_Y_ref;
        end
        
    catch ME
        % 如果发生错误，记录下来并返回标记值
        fprintf('  [FAIL] An FFT function call crashed: %s\n', ME.message);
        exec_time = NaN;
        result = [];
        error_norm = Inf;
    end
end