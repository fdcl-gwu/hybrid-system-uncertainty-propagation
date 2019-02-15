clear; close all;

% parameters
miu = -2;
theta = 0.5;
sigma = 0.2;
x0 = 1;
D = sigma^2/2;

% density function
f = @(t,x)sqrt(theta/(2*pi*D*(1-exp(-2*theta*t))))*exp(-theta/(2*D)*...
    (x-miu-(x0-miu)*exp(-theta*t)).^2/(1-exp(-2*theta*t)));

% grid
Lx = 10;
nx = 100;
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

% fft of f(x) = x
y_x = fftshift(fft(x))/nx;
y_x = y_x.*shift;

% reconstruct density
freq = ((0:nx-1)-floor(nx/2))/Lx;
f_fft = @(x,y)real(sum(y.'.*exp(1i*freq*2*pi.*x),2));

%% convolve first, then take derivation
y1(:,1) = ytrue(:,1);
A = zeros(nx,nx);
for n = 1:nx
    N = n-1-floor(nx/2);
    for k = 1:nx
        K = k-1-floor(nx/2);
        NMinusK = wrapDFT(N-K,nx);
        nMinusk = NMinusK+1+floor(nx/2);
        
        if n == 1 && mod(nx,2) == 0
            if K == 0
                A(n,nMinusk) = -2*sigma^2*pi^2*N^2/Lx^2;
            else
                A(n,nMinusk) = 0;
            end
        else
            if K == 0
                A(n,nMinusk) = -theta*miu*2*pi*1i*N/Lx + theta*2*pi*1i*N/Lx*y_x(k)...
                    -2*sigma^2*pi^2*N^2/Lx^2;
            else
                A(n,nMinusk) = theta*2*pi*1i*N/Lx*y_x(k);
            end
        end
    end
end

for i = 2:nt
    y1(:,i) = expm(A*(i-1)*Lt/nt)*y1(:,1);
end

%% take derivation first, then convolve
y2(:,1) = ytrue(:,1);
A = zeros(nx,nx);
for n = 1:nx
    N = n-1-floor(nx/2);
    for k = 1:nx
        K = k-1-floor(nx/2);
        NMinusK = wrapDFT(N-K,nx);
        nMinusk = NMinusK+1+floor(nx/2);
        
        if N == K
            if n == 1 && mod(nx,2) == 0
                part1 = theta-2*sigma^2*pi^2*N^2/Lx^2;
            else
                part1 = theta-theta*miu*2*pi*1i*N/Lx-2*sigma^2*pi^2*N^2/Lx^2;
            end
        else
            part1 = 0;
        end
        
        if k == 1 && mod(nx,2) == 0
            part2 = 0;
        else
            part2 = theta*2*pi*1i*K/Lx*y_x(nMinusk);
        end
        
        A(n,k) = part1+part2;
    end
end

for i = 2:nt
    y2(:,i) = expm(A*(i-1)*Lt/nt)*y2(:,1);
end

%% numerical integration of 

%% plot
figure; hold on;
for i = 1:nt
    plot3(x,ones(nx,1)*t(i),fx(:,i),'b');
    plot3(x,ones(nx,1)*t(i),f_fft(x,ytrue(:,i)),'r');
    plot3(x,ones(nx,1)*t(i),f_fft(x,y1(:,i)),'g');
end

figure; hold on;
for i = 1:nt
    plot3(x,ones(nx,1)*t(i),fx(:,i),'b');
    plot3(x,ones(nx,1)*t(i),f_fft(x,ytrue(:,i)),'r');
    plot3(x,ones(nx,1)*t(i),f_fft(x,y2(:,i)),'g');
end

