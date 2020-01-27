function [ fx, y ] = cont( mode )

close all;
addpath('..','..\..\lib');

% parameters
Vh = -10;
Vc = 10;
b = 15;

% grid
N1 = 120; N2 = 120;
L1 = 600; L2 = 600;
x1 = linspace(-L1/2,L1/2-L1/N1,N1);
x2 = linspace(-L2/2,L2/2-L2/N2,N2);
Nt = 401;
Lt = 40;
t = linspace(0,Lt,Nt);

% initial conditions
x1_0 = 180; x2_0 = 0;
sigma1_0 = 10; sigma2_0 = 10;

% initial distribution
fx = zeros(N1,N2,Nt);
fx(:,:,1) = 1/(sqrt(2*pi)*sigma1_0)*exp(-0.5*(reshape(x1,[],1)-x1_0).^2/sigma1_0^2) .*...
    (1/sqrt(2*pi)*sigma2_0*exp(-0.5*(reshape(x2,1,[])-x2_0).^2/sigma2_0^2));

% Fourier transform of initial distribution
y = zeros(N1,N2,Nt);
shift1 = reshape((-1).^((0:N1-1)-N1/2),[],1);
shift2 = reshape((-1).^((0:N2-1)-N2/2),1,[]);
y(:,:,1) = fftshift(fft2(fx(:,:,1)))/N1/N2;
y(:,:,1) = y(:,:,1).*shift1.*shift2;

% coefficients of approximated system
A = zeros(N1,N2);
c1 = [0,ones(1,N1-1)];
c2 = [0,ones(1,N2-1)];
for n1 = 1:N1
    for n2 = 1:N2
        part1 = -2*pi*1i*(n1-1-N1/2)*c1(n1)*Vh/L1;
        part2 = -2*pi*1i*(n2-1-N2/2)*c2(n2)*(mode==2)*Vc/L2;
        part3 = -0.5*b^2*(4*pi^2*(n1-1-N1/2)^2/L1^2 + 4*pi^2*(n2-1-N2/2)^2/L2^2);
        A(n1,n2) = part1+part2+part3;
    end
end

% exponential of the coefficient matrix
expA = zeros(N1,N2);
for n1 = 1:N1
    for n2 = 1:N2
        expA(n1,n2) = exp(A(n1,n2)*Lt/(Nt-1));
    end
end

% propagation
for nt = 2:Nt
    for n1 = 1:N1
        for n2 = 1:N2
            y(n1,n2,nt) = expA(n1,n2)*y(n1,n2,nt-1);
        end
    end
end

% Inverse Fourier transform
for nt = 1:Nt
    fx(:,:,nt) = ifft2(ifftshift(y(:,:,nt)./shift1./shift2*N1*N2),'symmetric');
end

for nt = 1:10:Nt
    figure;
    surf(x1,x2,fx(:,:,nt));
    view([0,0,1]);
end

rmpath('..','..\..\lib');

end

