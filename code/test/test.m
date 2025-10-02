clear;
clc;
close all;

fprintf('Benchmark is running silently. Please wait for the final results table...\n');
fprintf('This may take several minutes depending on the test range.\n\n');

% --- 1. 设置测试参数 ---
min_pow = 4;
max_pow = 10;
N_values = 4.^(min_pow:max_pow);

func_handles = {
    @fftNewPadded, @fft_iterative, @fft_radix4, ...
    @fft_radix4_iterative, @stockham_fft, @fftdivpad
};

func_names = {
    'R2-Recursive', 'R2-Iterative', 'R4-Recursive', ...
    'R4-Iterative', 'Stockham', 'Split-Radix'
};

% --- 2. 性能与正确性综合基准测试 (无打印输出) ---
num_algorithms = length(func_handles);
num_N_values = length(N_values);

results_time = zeros(num_N_values, num_algorithms);
results_error = zeros(num_N_values, num_algorithms);

for i = 1:num_N_values
    N = N_values(i);
    X = randn(1, N) + 1i * randn(1, N);
    Y_ref = fft(X);
    norm_Y_ref = norm(Y_ref);
    % 修改后的代码
    for j = 1:num_algorithms
        current_func_name = func_names{j};
        
        % --- 区分对待Stockham和其他函数 ---
        if strcmp(current_func_name, 'Stockham')
            % 如果是Stockham，使用列向量
            X_input = X.'; % 转置为列向量
            Y_test = func_handles{j}(X_input);
            Y_test = Y_test.'; % 转置回来以便比较
            f_to_test = @() func_handles{j}(X_input);
        else
            % 其他函数，使用行向量
            X_input = X;
            Y_test = func_handles{j}(X_input);
            f_to_test = @() func_handles{j}(X_input);
        end
    
        % --- 性能测试 (使用 f_to_test) ---
        try
            results_time(i, j) = timeit(f_to_test);
        catch
            results_time(i, j) = NaN;
        end
        
        % --- 正确性测试 (使用 Y_test) ---
        try
            if ~isequal(size(Y_test), size(Y_ref))
                error('Output dimension mismatch.');
            end
            results_error(i, j) = norm(Y_test - Y_ref) / norm_Y_ref;
        catch
            results_error(i, j) = Inf;
        end
    end
end

% --- 3. 格式化输出最终的单一表格 ---
% 合并时间和误差到一个cell数组中以便于格式化
num_cols_per_algo = 2; % Time, Error
table_data = cell(num_N_values, 1 + num_algorithms * num_cols_per_algo);

for i = 1:num_N_values
    table_data{i, 1} = N_values(i);
    for j = 1:num_algorithms
        col_idx = 1 + (j-1)*num_cols_per_algo;
        table_data{i, col_idx + 1} = results_time(i, j);
        table_data{i, col_idx + 2} = results_error(i, j);
    end
end

% 打印表头
fprintf('=================================== FFT Performance & Correctness Summary ===================================\n');
fprintf('%-10s', 'N');
for i = 1:num_algorithms
    fprintf('| %-20s ', func_names{i});
end
fprintf('|\n');

fprintf('%-10s', '');
for i = 1:num_algorithms
    fprintf('| %-10s %-10s', 'Time (s)', 'Error');
end
fprintf('|\n');

fprintf(repmat('=', 1, 10 + num_algorithms * 23));
fprintf('\n');

% 打印数据行
for i = 1:num_N_values
    fprintf('%-10d', table_data{i,1});
    for j = 2:size(table_data, 2)
        if mod(j, 2) == 0 % Time column
            fprintf('| %-10.6f ', table_data{i,j});
        else % Error column
            fprintf('%-10.2e', table_data{i,j});
        end
    end
    fprintf('|\n');
end
fprintf(repmat('=', 1, 10 + num_algorithms * 23));
fprintf('\n');
fprintf('\nBenchmark finished. Data table is ready to be copied.\n');