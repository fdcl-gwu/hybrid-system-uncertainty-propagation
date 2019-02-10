clear; close all;

% parameters
g = 9.8;
niu = 5;
sigma = 0.5;
x0 = [8;0];
sigma0 = [0.5^2,0;0,0.1^2];

% grid
n1 = 100; n2 = 100;
L1 = 20; L2 = 10;
x1 = linspace(-L1/2,L1/2-L1/n1,n1);
x2 = linspace(-L2/2,L2/2-L2/n2,n2);
nt = 21;
Lt = 5;
t = linspace(0,Lt,nt);

% initial density function
f0 = @(x)1/(2*pi)/sqrt(det(sigma0))*exp(-1/2*(x-x0)'*sigma0^-1*(x-x0));
fx = zeros(n1,n2,nt);
for i = 1:n1
    for j = 1:n2
        fx(i,j,1) = f0([x1(i);x2(j)]);
    end
end

% initial fft
y = zeros(n1,n2,nt);
y(:,:,1) = fftshift(fft2(fx(:,:,1)))/n1/n2;

% fft for f(x)=x and f(x)=x^2
y_x = fftshift(fft(x2))/n2;
y_xsqr = fftshift(fft(x2.^2))/n2;

% propagation
y0 = y(:,:,1);
parfor m = 1:n1
    A = zeros(n2,n2);
    M = m-1-floor(n1/2);
    for i = 2:nt
        for n = 1:n2
            N = n-1-floor(n2/2);
            for k = 1:n2
                K = k-1-floor(n2/2);
                NMinusK = wrapDFT(N-K,n2);
                nMinusk = NMinusK+1+floor(n2/2);
                
                if m == 1 && mod(n1,2) == 0
                    part1 = 0;
                else
                    part1 = -2*pi*1i*M/L1*y_x(k);
                end
                
                if n == 1 && mod(n2,2) == 0
                    part2 = 0;
                else
                    if K == 0
                        part2 = 2*pi*1i*N/L2*g + 2*pi*1i*N/L2*niu*y_x(k);
                    else
                        part2 = 2*pi*1i*N/L2*niu*y_x(k);
                    end
                end
                
                part3 = -2*sigma^2*pi^2*N^2/L2^2*y_xsqr(k);
                
                A(n,nMinusk) = part1 + part2 + part3;
            end
        end
        y(m,:,i) = (expm(A*(i-1)*Lt/(nt-1))*y0(m,:).').';
    end
end

% reconstruct density
fraq1 = ((0:n1-1)'-floor(n1/2))/L1;
fraq2 = ((0:n2-1)-floor(n2/2))/L2;

f_fft = @(x,y)sum(sum(y.*exp(1i*fraq1*2*pi*x(1)).*exp(1i*fraq2*2*pi*x(2))));
parfor k = 2:nt
    for i = 1:n1
        for j = 1:n2
            fx(i,j,k) = real(f_fft([x1(i)+L1/2;x2(j)+L2/2],y(:,:,k)));
        end
    end
end

% plot
for i = 1:nt
figure;
surf(x2,x1,fx(:,:,i));
view([0,0,1]);
end

