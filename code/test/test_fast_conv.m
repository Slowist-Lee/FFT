% --- test_fast_conv.m ---
clear; clc;

% 1. 创建两个随机信号
f = randn(1, 1000);
g = randn(1, 800);

% 2. 设定混合算法的切换点 (使用你找到的最佳值)
best_crossover_n = 416; 

% 3. 使用我们的快速卷积函数
tic;
y_fast = fast_conv(f, g, best_crossover_n);
toc;

% 4. 使用MATLAB内置的卷积函数作为对比
tic;
y_ref = conv(f, g);
toc;

% 5. 验证结果的正确性
error = norm(y_fast - y_ref) / norm(y_ref);
fprintf('Relative error between fast_conv and built-in conv: %e\n', error);

% 6. 绘制结果对比图
figure;
plot(y_fast, 'b-', 'LineWidth', 1.5);
hold on;
plot(y_ref, 'r--', 'LineWidth', 1.5);
hold off;
grid on;
title('Comparison of fastconv and built-in conv');
legend('Our fastconv Result', 'MATLAB conv() Result');
xlabel('Sample Index');
ylabel('Amplitude');