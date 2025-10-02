% =========================================================================
% 寻找 Stockham 与 Split-Radix 混合算法的最佳切换点 (CROSSOVER_N)
% =========================================================================

clear;
clc;
close all;

% --- 1. 全局设置 ---

% 定义用于评估性能的FFT规模范围
% 我们将测试混合算法在处理这一系列不同规模N时的平均性能
FFT_SIZES_TO_TEST = 2.^(10:20); % 测试从 N=1024 到 N=1048576 的性能

fprintf('Performance will be evaluated over FFT sizes from N=%d to N=%d.\n', ...
    min(FFT_SIZES_TO_TEST), max(FFT_SIZES_TO_TEST));
fprintf('This may take a few minutes...\n\n');


% 阶段一：粗略搜索
fprintf('--- Stage 1: Coarse Search ---\n');

coarse_crossover_points = 2.^(8:16); 
avg_times_coarse = zeros(size(coarse_crossover_points));

for i = 1:length(coarse_crossover_points)
    crossover_n = coarse_crossover_points(i);
    fprintf('  Testing coarse crossover point N = %d...\n', crossover_n);
    
    total_time = 0;
    for n_fft = FFT_SIZES_TO_TEST
        X = randn(1, n_fft) + 1i * randn(1, n_fft);
        f_hybrid = @() hybrid_fft(X, crossover_n);
        total_time = total_time + timeit(f_hybrid);
    end
    avg_times_coarse(i) = total_time / length(FFT_SIZES_TO_TEST);
end

figure;
semilogx(coarse_crossover_points, avg_times_coarse, 'b-o', 'LineWidth', 1.5);
grid on;
title('Stage 1: Coarse Search for Best Crossover Point');
xlabel('Crossover Point N (Log Scale)');
ylabel('Average Execution Time (s)');
xticks(coarse_crossover_points);
xtickangle(45);

[min_time_coarse, idx_coarse] = min(avg_times_coarse);
best_coarse_n = coarse_crossover_points(idx_coarse);
fprintf('\nCoarse search complete. Best performance found around N = %d.\n\n', best_coarse_n);


%阶段二：精细搜索
fprintf('--- Stage 2: Fine Search ---\n');

fine_search_range_start = max(2, best_coarse_n / 2);
fine_search_range_end = best_coarse_n * 2;
fine_search_step = 2^(log2(best_coarse_n)-3);
if fine_search_step < 32; fine_search_step=32; end

fine_crossover_points = fine_search_range_start : fine_search_step : fine_search_range_end;

fprintf('Performing fine search around N=%d, from %d to %d with step %d.\n', ...
    best_coarse_n, fine_search_range_start, fine_search_range_end, fine_search_step);

avg_times_fine = zeros(size(fine_crossover_points));

for i = 1:length(fine_crossover_points)
    crossover_n = fine_crossover_points(i);
    fprintf('  Testing fine crossover point N = %d...\n', crossover_n);
    
    total_time = 0;
    for n_fft = FFT_SIZES_TO_TEST
        X = randn(1, n_fft) + 1i * randn(1, n_fft);
        f_hybrid = @() hybrid_fft(X, crossover_n);
        total_time = total_time + timeit(f_hybrid);
    end
    avg_times_fine(i) = total_time / length(FFT_SIZES_TO_TEST);
end

% 绘制结果
figure;
plot(fine_crossover_points, avg_times_fine, 'r-s', 'LineWidth', 1.5);
grid on;
title('Stage 2: Fine Search for Best Crossover Point');
xlabel('Crossover Point N');
ylabel('Average Execution Time (s)');

[min_time_fine, idx_fine] = min(avg_times_fine);
best_fine_n = fine_crossover_points(idx_fine);

hold on;
plot(best_fine_n, min_time_fine, 'gp', 'MarkerSize', 15, 'MarkerFaceColor', 'g');
legend('Performance', sprintf('Best N = %d', best_fine_n), 'Location', 'best');
hold off;

fprintf('\n=========================================================\n');
fprintf('                OPTIMAL CROSSOVER POINT FOUND\n');
fprintf('=========================================================\n');
fprintf('The best crossover point N for your system is: %d\n', best_fine_n);
fprintf('=========================================================\n');


