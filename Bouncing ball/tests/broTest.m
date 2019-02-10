clear; close all;

% parameters
miu = -0.5;
sigma = 1;

% density function
f = @(t,x)1/sigma/sqrt(2*pi*t)*exp(-(x-miu*t).^2/(2*sigma^2*t));

% grid
Lx = 20;
nx = 500; N = ceil(nx/2)-1;
x = linspace(-Lx/2,Lx/2,nx)';
Lt = 5;
nt = 20;
t = linspace(Lt/nt,Lt,nt);

% true density values
fx = zeros(nx,nt);
for i = 1:nt
    fx(:,i) = f(t(i),x);
end

% true fft
ytrue = zeros(nx,nt);
for i = 1:nt
    ytrue(:,i) = fftshift(fft(fx(:,i)))/nx;
end
if mod(nx,2) == 0
    ytrue(1,:) = [];
end

% propagation
y(:,1) = ytrue(:,1);
A = zeros(2*N+1,2*N+1);
for i = 1:2*N+1
    I = i-1-N;
    A(i,i) = -miu*2*pi*1i*I/Lx-2*sigma^2*pi^2*I^2/Lx^2;
end

for i = 2:nt
    y(:,i) = expm(A*(i-1)*Lt/nt)*y(:,1);
end

% reconstruct density
freq = (-N:N)/Lx*2*N/(2*N+1);
f_fft = @(x,y)real(sum(y.'.*exp(1i*freq*2*pi.*x),2));

% plot
figure; hold on;
for i = 1:nt
    plot3(x,ones(nx,1)*t(i),fx(:,i),'b');
    plot3(x,ones(nx,1)*t(i),f_fft(x+Lx/2,ytrue(:,i)),'r');
    plot3(x,ones(nx,1)*t(i),f_fft(x+Lx/2,y(:,i)),'g');
end


