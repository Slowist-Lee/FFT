function y = fast_conv(f, g, crossover_n)
    Nf = length(f);
    Ng = length(g);
    %计算线性卷积结果的长度
    L = Nf + Ng - 1;
    N_fft = 2^nextpow2(L);
    fp = [f, zeros(1, N_fft - Nf)];
    gp = [g, zeros(1, N_fft - Ng)];
    %用hybrid_fft正向FFT
    Fp = hybrid_fft(fp, crossover_n);
    Gp = hybrid_fft(gp, crossover_n);
    Yp = Fp .* Gp;
    %傅里叶逆变换
    Y_conj_fft = hybrid_fft(conj(Yp), crossover_n);
    yp = (1/N_fft) * conj(Y_conj_fft);
    %最终的线性卷积结果, 取实部并截取前 L 个点
    y = real(yp(1:L));
end