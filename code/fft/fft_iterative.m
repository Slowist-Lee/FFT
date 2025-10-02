function Y = fft_iterative(X)
    %和之前一样先补零
    N=length(X);
    p=nextpow2(N);
    L=2^p;
    if N==L
        X_padded=X;
    else
        X_padded=[X,zeros(1,L-N)]; %在X后补L-N个零
    end
    
    X_reversed=zeros(1,L);
    for i=1:L
        X_reversed(reverse(i-1,p)+1)=X_padded(i);
    end
    Y=X_reversed;
    %m代表FFT的规模，逐渐变大
    m=2;
    while m<=L
        h=m/2;
        %旋转因子
        dw=exp(-1j*2*pi/m);
        %对整个序列以m个为一组进行蝶形计算
        for k=1:m:L
            w=1;
            for j=0:h-1
                i1=k+j;
                i2=k+j+h;
                %一个小组里的蝶形计算
                temp = Y(i1);
                term = w * Y(i2);
                Y(i1) = temp + term;
                Y(i2) = temp - term;
                w=w*dw; %更新旋转因子
            end
        end
        m=m*2;
    end
end


%实现位反转
function r = reverse(n,p)
    r=0;
    for i=1:p
        l=mod(n,2);
        r=r*2+l;
        n=floor(n/2);
    end
end
