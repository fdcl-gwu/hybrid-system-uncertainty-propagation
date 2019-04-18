function [fx] = cont()
close all;
addpath('..','..\..\lib');

p = getParameter(1);
% parameters
g = p.g;
niu = p.niu;
sigma = p.sigma;
x0 = p.x0;
sigma0 = p.sigma0;

% grid
N1 = p.N1; N2 = p.N2;
L1 = p.L1; L2 = p.L2;
x1 = linspace(-L1/2,L1/2-L1/N1,N1);
x2 = linspace(-L2/2,L2/2-L2/N2,N2);
Nt = 41;
Lt = 1;
t = linspace(0,Lt,Nt);

% initial density function
f0 = @(x)1/(2*pi)/sqrt(det(sigma0))*exp(-1/2*(x-x0)'*sigma0^-1*(x-x0));
fx = zeros(N1,N2,Nt);
for n1 = 1:N1
    for n2 = 1:N2
        fx(n1,n2,1) = f0([x1(n1);x2(n2)]);
    end
end

% initial fft
shift1 = (-1).^((0:N1-1)-floor(N1/2)).';
shift2 = (-1).^((0:N2-1)-floor(N2/2));
y = zeros(N1,N2,Nt);
y(:,:,1) = fftshift(fft2(fx(:,:,1)))/N1/N2;
y(:,:,1) = shift1.*shift2.*y(:,:,1);

% fft for f(x)=x and f(x)=x^2
y_x = fftshift(fft(x2))/N2;
y_x = shift2.*y_x;
y_xsqr = fftshift(fft(x2.*abs(x2)))/N2;
y_xsqr = shift2.*y_xsqr;
y_xfour = fftshift(fft(x2.^4))/N2;
y_xfour = shift2.*y_xfour;

% coefficients of approximated systems
A = zeros(N2,N2,N1);
for n1 = 1:N1
    for n2 = 1:N2
        for m2 = 1:N2
            n2Minusm2 = wrapDFT(n2-m2,N2)+1+N2/2;
            
            if n1 == 1
                part1 = 0;
            else
                part1 = -2*pi*1i*(n1-1-N1/2)/L1*y_x(m2);
            end
            
            if n2 == 2
                part2 = 0;
            else
                if m2-1-N2/2 == 0
                    part2 = 2*pi*1i*(n2-1-N2/2)/L2*g + 2*pi*1i*(n2-1-N2/2)/L2*niu*y_xsqr(m2);
                else
                    part2 = 2*pi*1i*(n2-1-N2/2)/L2*niu*y_xsqr(m2);
                end
            end
            
            part3 = -2*sigma^2*pi^2*(n2-1-N2/2)^2/L2^2*y_xfour(m2);
            
            A(n2,n2Minusm2,n1) = part1 + part2 + part3;
        end
    end
end

expA = zeros(N2,N2,N1);
for n1 = 1:N1
    expA(:,:,n1) = expm(A(:,:,n1)*Lt/(Nt-1));
end
clear A;

% propagation
for nt = 2:Nt
    for n1 = 1:N1
        y(n1,:,nt) = (expA(:,:,n1)*y(n1,:,nt-1).').';
    end
end

% reconstruct density
parfor nt = 2:Nt
    fx(:,:,nt) = ifftn(ifftshift(y(:,:,nt)./shift1./shift2*N1*n2),'symmetric');
end

% plot
for n1 = 1:Nt
    figure;
    surf(x2,x1,fx(:,:,n1));
    view([0,0,1]);
end

rmpath('..','..\..\lib');

end

