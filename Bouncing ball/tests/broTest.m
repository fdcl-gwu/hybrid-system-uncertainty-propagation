clear; close all;

% parameters
miu = -0.5;
sigma = 1;

% density function
f = @(t,x)1/sigma/sqrt(2*pi*t)*exp(-(x-miu*t).^2/(2*sigma^2*t));

% grid
Lx = 20;
nx = 500;
x = linspace(-Lx/2,Lx/2-Lx/nx,nx)';
Lt = 5;
nt = 20;
t = linspace(Lt/nt,Lt,nt);

% true density values
fx = zeros(nx,nt);
for i = 1:nt
    fx(:,i) = f(t(i),x);
end

% true fft
shift = (-1).^((0:nx-1)-floor(nx/2)).';
ytrue = zeros(nx,nt);
for i = 1:nt
    ytrue(:,i) = fftshift(fft(fx(:,i)))/nx;
end
ytrue = ytrue.*shift;

% propagation
y(:,1) = ytrue(:,1);
A = zeros(nx,nx);
for n = 1:nx
    N = n-1-floor(nx/2);
    if n == 1 && mod(nx,2) == 0
        A(n,n) = -2*sigma^2*pi^2*N^2/Lx^2;
    else
        A(n,n) = -miu*2*pi*1i*N/Lx-2*sigma^2*pi^2*N^2/Lx^2;
    end
end

for i = 2:nt
    y(:,i) = expm(A*(i-1)*Lt/nt)*y(:,1);
end

% reconstruct density
freq = ((0:nx-1)-floor(nx/2))/Lx;
f_fft = @(x,y)real(sum(y.'.*exp(1i*freq*2*pi.*x),2));

% plot
figure; hold on;
for i = 1:nt
    plot3(x,ones(nx,1)*t(i),fx(:,i),'b');
    plot3(x,ones(nx,1)*t(i),f_fft(x,ytrue(:,i)),'r');
    plot3(x,ones(nx,1)*t(i),f_fft(x,y(:,i)),'g');
end


