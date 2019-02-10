clear;
close all;

% geometric Brownian motion parameters
miu = 0;
sigma = 0.5;
x0 = 3;

% GBM density function
f = @(t,x)1/sqrt(2*pi)./(x*sigma*sqrt(t)).*exp(-(log(x)-log(x0)-t*(miu-1/2*sigma^2)).^2/(2*sigma^2*t));

% coordinates
Lx = 5;
nx = 101; N = ceil(nx/2)-1;
x = linspace(0,Lx,nx)';
Lt = 5;
nt = 20;
t = linspace(Lt/nt,Lt,nt);

% true density values
fx = zeros(nx,nt);
for i = 1:nt
    fx(1,i) = 0;
    fx(2:nx,i) = f(t(i),x(2:nx));
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
    for j = 1:2*N+1
        J = j-1-N;
        K = wrapToN(I-J,N);
        k = K+1+N;
        if J == 0
            A(i,k) = -sigma^2*pi^2*I^2/3*2;
        else
            A(i,k) = miu*I/J-sigma^2*I^2*(1/J^2+1i*pi/J);
        end
    end
end

for i = 2:nt
    y(:,i) = expm(A*(i-1)*Lt/nt)*y(:,1);
end

% reconstruct
fraq = (-N:N)/Lx*2*N/(2*N+1);
f_fft = @(x,y)real(sum(y.'.*exp(1i*fraq*2*pi.*x),2));

% plot
figure; hold on;
for i = 1:nt
    plot3(x,ones(nx,1)*t(i),fx(:,i),'b');
    plot3(x,ones(nx,1)*t(i),f_fft(x,ytrue(:,i)),'r');
    plot3(x,ones(nx,1)*t(i),f_fft(x,y(:,i)),'g');
end

