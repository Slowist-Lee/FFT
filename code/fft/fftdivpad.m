function Y = fftdivpad(X)
    N=length(X);
    p=nextpow2(N);
    L=2^p;
    if N==L
        X_padded=X;
    else
        X_padded=[X,zeros(1,L-N)]; %在X后补L-N个零
    end
    Y=zeros(1,L);
    Y=fftDiv(X_padded);
end

function Y = fftDiv(X)
    if length(X) == 2
        Y = zeros(1,2);
        Y(1) = X(1)+X(2);
        Y(2) = X(1)-X(2);
    elseif length(X)==1
        Y=zeros(1,1);
        Y(1)=X(1);
    else
        N=length(X);
        X1 = X([1:2:N]);%假设N=16,那么第一句取了x[1,3,5,7,9,...15],0开始就是0,2,4,...14
        X2 = X([2:4:N]);%取1,5,9,13
        X3 = X([4:4:N]);%3,7,11,15
        Y1 = fftDiv(X1);
        Y2 = fftDiv(X2);
        Y3 = fftDiv(X3);%递归。
        %下面这一段讲的就是合并
        Y = zeros(1,N);

        for k=1:N/4
            T1=Y2(k)*exp(-1j*2*pi*(k-1)/N);
            T2=Y3(k)*exp(-1j*2*pi*3*(k-1)/N);
            Y(k)=Y1(k)+T1+T2;
            Y(k+N/2)=Y1(k)-T1-T2;
            Y(k+N/4)=Y1(k+N/4)-1j*(T1-T2);
            Y(k+3*N/4)=Y1(k+N/4)+1j*(T1-T2);
        end
    end
end
