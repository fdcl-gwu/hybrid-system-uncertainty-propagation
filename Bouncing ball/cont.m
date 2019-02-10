clear; close all;

% parameters
g = 9.8;
niu = 5;
sigma = 0.5;
x0 = [8;0];
sigma0 = [0.5^2,0;0,0.1^2];

% grid
n1 = 100; n2 = 100;
N1 = ceil(n1/2)-1; N2 = ceil(n2/2)-1;
L1 = 20; L2 = 10;
x1 = linspace(-L1/2,L1/2,n1);
x2 = linspace(-L2/2,L2/2,n2);
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
if mod(n1,2) == 0
    y(1,:,:) = [];
end
if mod(n2,2) == 0
    y(:,1,:) = [];
end

% propagation
y0 = y(:,:,1);
for m = 1:2*N1+1
    A = zeros(2*N2+1,2*N2+1);
    M = m-N1-1;
    for i = 2:nt
        for n = 1:2*N2+1
            N = n-N2-1;
            for k = 1:2*N2+1
                K = k-N2-1;
                NMinusK = wrapToN(N-K,N2);
                nMinusk = NMinusK+N2+1;
                if K == 0
                    A(n,n) = g*2*pi*1i*N/L2-sigma^2*pi^2*N^2/6;
                else
                    A(n,nMinusk) = M*L2/L1/K-niu*N/K-sigma^2*N^2/K^2;
                end
            end
        end
        y(m,:,i) = (expm(A*(i-1)*Lt/(nt-1))*y0(m,:).').';
    end
end

% reconstruct density
fraq1 = (-N1:N1)'/L1*2*N1/(2*N1+1);
fraq2 = (-N2:N2)/L2*2*N2/(2*N2+1);

f_fft = @(x,y)sum(sum(y.*exp(1i*fraq1*2*pi*x(1)).*exp(1i*fraq2*2*pi*x(2))));
for k = 2:nt
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

