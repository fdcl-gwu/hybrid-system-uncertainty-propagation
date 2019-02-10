clear;
close all;

% Gaussian distribution parameters
miu = 2;
sigma = 1;

% density function
f = @(x)1/sqrt(2*pi)/sigma*exp(-(x-miu).^2/2/sigma^2);

% grid
L = 10;
n = 1000;
x = linspace(-L/2,L/2-L/n,n)';

% true density value
fx = f(x);

% fft
y = fftshift(fft(fx))/n;

% reconstruct density function
freq = ((0:n-1)-floor(n/2))/L;
f_fft = @(x)real(sum(y.'.*exp(1i*freq*2*pi.*x),2));

% plot
n = 1001;
x = linspace(-L/2,L/2-L/n,n)';

figure; hold on;
plot(x,f(x));
plot(x,f_fft(x+L/2));

