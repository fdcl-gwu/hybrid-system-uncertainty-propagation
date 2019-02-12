clear;
close all;

% geometric Brownian motion parameters
miu = -0.2;
sigma = 0.5;
x0 = 2;

% GBM density function
f = @(t,x)1/sqrt(2*pi)./(x*sigma*sqrt(t)).*exp(-(log(x)-log(x0)-t*(miu-1/2*sigma^2)).^2/(2*sigma^2*t));

% coordinates
Lx = 10;
nx = 500;
x = linspace(-Lx/2,Lx/2-Lx/nx,nx)';
Lt = 5;
nt = 20;
t = linspace(Lt/nt,Lt,nt);

% true density values
fx = zeros(nx,nt);
for i = 1:nt
    fx(x<=0) = 0;
    fx(x>0,i) = f(t(i),x(x>0));
end

% true fft
shift = (-1).^((0:nx-1)-floor(nx/2)).';
ytrue = zeros(nx,nt);
for i = 1:nt
    ytrue(:,i) = fftshift(fft(fx(:,i)))/nx;
end
ytrue = ytrue.*shift;

% fft of f(x)=x and f(x)=x^2;
y_x = fftshift(fft(x))/nx;
y_x = y_x.*shift;
y_xsqr = fftshift(fft(x.^2))/nx;
y_xsqr = y_xsqr.*shift;

% propagation
y(:,1) = ytrue(:,1);
A = zeros(nx,nx);
for n = 1:nx
    N = n-1-floor(nx/2);
    for k = 1:nx
        K = k-1-floor(nx/2);
        NMiunsK = wrapDFT(N-K,nx);
        nMinusk = NMiunsK+1+floor(nx/2);
        
        if n == 1 && mod(nx,2) == 0
            A(n,nMinusk) = -2*sigma^2*pi^2*N^2/Lx^2*y_xsqr(k);
        else
            A(n,nMinusk) = -miu*2*pi*1i*N/Lx*y_x(k) - 2*sigma^2*pi^2*N^2/Lx^2*y_xsqr(k);
        end
    end
end

for i = 2:nt
    y(:,i) = expm(A*(i-1)*Lt/nt)*y(:,1);
end

% reconstruct
fraq = ((0:nx-1)-floor(nx/2))/Lx;
f_fft = @(x,y)real(sum(y.'.*exp(1i*fraq*2*pi.*x),2));

% plot
figure; hold on;
for i = 1:nt
    plot3(x,ones(nx,1)*t(i),fx(:,i),'b');
    plot3(x,ones(nx,1)*t(i),f_fft(x,ytrue(:,i)),'r');
    plot3(x,ones(nx,1)*t(i),f_fft(x,y(:,i)),'g');
end

