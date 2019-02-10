clear; close all;

% parameters
miu = -2;
theta = 0.5;
sigma = 1;
x0 = 1;
D = sigma^2/2;

% density function
f = @(t,x)sqrt(theta/(2*pi*D*(1-exp(-2*theta*t))))*exp(-theta/(2*D)*...
    (x-miu-(x0-miu)*exp(-theta*t)).^2/(1-exp(-2*theta*t)));

% grid
Lx = 10;
nx = 101; N = ceil(nx/2)-1;
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
    for j = 1:2*N+1
        J = j-1-N;
        K = wrapToN(I-J,N);
        k = K+1+N;
        if J == 0
            A(i,k) = -2*sigma^2*pi^2*I^2/Lx^2-miu*theta*2*pi*1i*I/Lx;
        else
            A(i,k) = -theta*I/J;
        end
    end
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

