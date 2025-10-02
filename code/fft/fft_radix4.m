function Y = fft_radix4(X)
    N=length(X);
    %先补零
    N=length(X);
    p=ceil(log2(N)/2);
    L=4^p;
    if N==L
        X_padded=X;
    else
        X_padded=[X,zeros(1,L-N)]; %在X后补L-N个零
    end
    Y=zeros(1,L);
    Y = fft_radix4rec(X_padded);
end

function Y = fft_radix4rec(X)
    N = length(X);
    
    % 递归终止条件：4点FFT
    if N == 4
        Y = zeros(1, 4);
        Y(1) = X(1) + X(2) + X(3) + X(4);
        Y(2) = X(1) - 1j*X(2) - X(3) + 1j*X(4);
        Y(3) = X(1) - X(2) + X(3) - X(4);
        Y(4) = X(1) + 1j*X(2) - X(3) - 1j*X(4);
        return;
    end
    
    % 1. 分解成4个子序列
    X0 = X(1:4:N);   % X[0], X[4], X[8], ...
    X1 = X(2:4:N);   % X[1], X[5], X[9], ...
    X2 = X(3:4:N);   % X[2], X[6], X[10], ...
    X3 = X(4:4:N);   % X[3], X[7], X[11], ...
    
    % 2. 四个子序列递归进行FFT运算
    Y0 = fft_radix4rec(X0);  % b_k
    Y1 = fft_radix4rec(X1);  % c_k
    Y2 = fft_radix4rec(X2);  % d_k
    Y3 = fft_radix4rec(X3);  % e_k
    
    % 3. 蝶形合并
    Y = zeros(1, N);
    
    % 对每个k进行合并 (k从0开始，但MATLAB索引从1开始)
    for k = 0:(N/4-1)
        % 计算旋转因子
        W1 = exp(-1j * 2 * pi * k / N);      % W_N^k
        W2 = exp(-1j * 2 * pi * 2 * k / N);  % W_N^(2k)
        W3 = exp(-1j * 2 * pi * 3 * k / N);  % W_N^(3k)
        
        % 索引转换：k -> k+1 (MATLAB索引)
        idx = k + 1;
        
        % 按照推导的合并公式计算4个输出点
        Y(idx)= Y0(idx) + W1 * Y1(idx) + W2 * Y2(idx) + W3 * Y3(idx);
        Y(idx + N/4)  = Y0(idx) + W1 * (-1j) * Y1(idx) + W2 * (-1) * Y2(idx) + W3 * (1j) * Y3(idx);
        Y(idx + N/2)  = Y0(idx) + W1 * (-1) * Y1(idx) + W2 * (1) * Y2(idx) + W3 * (-1) * Y3(idx);
        Y(idx + 3*N/4) = Y0(idx) + W1 * (1j) * Y1(idx) + W2 * (-1) * Y2(idx) + W3 * (-1j) * Y3(idx);
    end
end
