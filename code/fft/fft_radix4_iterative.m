function Y = fft_radix4_iterative(X)
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
    X_reversed=zeros(1,L);
    for i=1:L
        X_reversed(reverse(i-1,p)+1)=X_padded(i);
    end
    Y=X_reversed;
    %m代表FFT的规模，逐渐变大
    m=4;
    while m<=L
        q=m/4;
        % 预计算当前阶段 m 所需的所有旋转因子
        W_m = zeros(1, q);
        for j = 0 : q-1
            W_m(j+1) = exp(-1j * 2 * pi * j / m);
        end
        
        % 对整个序列以m个为一组进行蝶形计算
        for k=1:m:L
            for j=0:q-1
                % 不再计算exp()，而是直接从表中查找！
                w1 = W_m(j+1); 
                w2 = w1*w1;
                w3 = w2*w1;

                x0=Y(k+j);
                x1=Y(k+j+q);
                x2=Y(k+j+2*q);
                x3=Y(k+j+3*q);
                
                T1=w1*x1;
                T2=w2*x2;
                T3=w3*x3;

                Y(k+j) = x0+T1+T2+T3;
                Y(k+j+q) = x0 - 1j*T1 - T2 + 1j*T3;
                Y(k+j+2*q) = x0 - T1 + T2 - T3;
                Y(k+j+3*q) = x0 + 1j*T1 - T2 - 1j*T3;
            end
        end
        m=m*4;
    end
end


%实现位反转
function r = reverse(n,p)
    r=0;
    for i=1:p
        l=mod(n,4);
        r=r*4+l;
        n=floor(n/4);
    end
end