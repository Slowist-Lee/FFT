% --- 混合FFT函数定义 ---
% 这个函数被上面的测试循环所调用
function Y = hybrid_fft(X, crossover_n)
    % HYBRID_FFT: 根据crossover_n在Stockham和Radix-4-iterative间切换
    % 假设：N > crossover_n 使用 Stockham (适合大规模)
    %      N <= crossover_n 使用 Split-Radix (适合小规模)
    % !!! 注意: 请根据你的实际结论调整这个逻辑 !!!
    % 如果你的结论是反过来的，请交换下面if/else里的函数调用

    N = length(X);

    if N > crossover_n
        % 大规模 -> 调用Stockham
        % Stockham函数期望列向量，需要转置
        Y = stockham_fft(X.'); 
        Y = Y.'; % 转置回来
    else
        % 小规模 -> 调用radix4
        Y = fft_radix4_iterative(X);
    end
end