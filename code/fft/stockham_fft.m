function y = stockham_fft(x)
    % --- 1. 初始化和输入检查 ---
    x = x(:); % 确保x是一个列向量
    N_orig = length(x); % 保存原始长度以备后用（如果需要）
    
    % 计算补零后的长度 L (必须是2的幂)
    p=nextpow2(N_orig);
    L = 2^p; 

    % 如果需要，进行补零
    if N_orig < L
        x_padded = [x; zeros(L - N_orig, 1, 'like', x)]; % 正确的垂直补零
    else
        x_padded = x;
    end
    x=x_padded;
    N=L;
    % --- 2. 设置 "乒乓" 缓冲区 ---
    buffers = zeros(N, 2, 'like', x);
    buffers(:, 1) = x; % 将输入数据复制到第一个缓冲区
    in_idx = 1;  % 当前输入缓冲区的列索引 (1 或 2)
    out_idx = 2; % 当前输出缓冲区的列索引 (2 或 1)
    % --- 3. 迭代执行 FFT 阶段 ---
    % 总共有 log2(N) 个阶段
    for stage = 0:(p - 1)
        s = 2^stage;
        n = N / s;
        m = n / 2;
        % 对所有子问题进行蝶形运算
        for p = 0:(m - 1)
            % 计算旋转因子 W_n^p = exp(-j*2*pi*p/n)
            wp = exp(-1i * 2 * pi * p / n);
            for q = 0:(s - 1)
                idx_a = q + s*p + 1;
                idx_b = q + s*(p+m) + 1;
                
                a = buffers(idx_a, in_idx);
                b = buffers(idx_b, in_idx);
                
                % 蝶形运算 (DIF - Decimation In Frequency)
                % out_buf[even] = a + b
                % out_buf[odd]  = (a - b) * W
                idx_out_even = q + s*(2*p) + 1;
                idx_out_odd  = q + s*(2*p+1) + 1;

                buffers(idx_out_even, out_idx) = a + b;
                buffers(idx_out_odd, out_idx)  = (a - b) * wp;
            end
        end
        
        % "乒乓": 交换输入和输出缓冲区的角色，为下一阶段做准备
        temp_idx = in_idx;
        in_idx = out_idx;
        out_idx = temp_idx;
    end
    % --- 4. 返回最终结果 ---
    % 经过 log2(N) 个阶段后, 最终结果位于 'in_idx' 指向的缓冲区中
    y = buffers(:, in_idx);
end