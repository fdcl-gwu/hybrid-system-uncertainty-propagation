function [ fx, y ] = cont( mode )

close all;
addpath('..\lib');

% parameters
v = 0.5;
u = [0,0.5,-0.5];
sigma = 0.2;
if ~exist('mode','var') || isempty('mode')
    mode = 1;
end

% grid
N1 = 50; N2 = 50;
L1 = 6; L2 = 6;
x1 = linspace(-L1/2,L1/2-L1/N1,N1);
x2 = linspace(-L2/2,L2/2-L2/N2,N2);
N3 = 50;
x3 = linspace(-pi,pi-2*pi/N3,N3);
Nt = 41;
Lt = 4;
t = linspace(0,Lt,Nt);

% initial conditions
x1_0 = 0; x2_0 = -2;
sigma1_0 = 0.2; sigma2_0 = 0.2;
x3_0 = pi/2;
k_0 = 20;

% initial distribution
fx = zeros(N1,N2,N3,Nt);
fx(:,:,:,1) = 1/(sqrt(2*pi)*sigma1_0)*exp(-0.5*(reshape(x1,[],1,1)-x1_0).^2/sigma1_0^2) .* ...
    (1/(sqrt(2*pi)*sigma2_0)*exp(-0.5*(reshape(x2,1,[],1)-x2_0).^2/sigma2_0^2)) .* ...
    (1/(2*pi*besseli(0,k_0))*exp(k_0*cos(reshape(x2,1,1,[])-x3_0)));

% Fourier transform of initial distribution
y = zeros(N1,N2,N3,Nt);
shift1 = reshape((-1).^((0:N1-1)-N1/2),[],1,1);
shift2 = reshape((-1).^((0:N2-1)-N2/2),1,[],1);
shift3 = reshape((-1).^((0:N3-1)-N3/2),1,1,[]);
y(:,:,:,1) = fftshift(fftn(fx(:,:,:,1)))/N1/N2/N3;
y(:,:,:,1) = y(:,:,:,1).*shift1.*shift2.*shift3;

% Fourier transfor for sin(theta) and cos(theta)
ysin = fftshift(fft(sin(x3)))/N3;
ysin = ysin.*reshape(shift3,1,[]);
ycos = fftshift(fft(cos(x3)))/N3;
ycos = ycos.*reshape(shift3,1,[]);

% coefficients of approximated system
A = zeros(N3,N3,N1,N2);
c1 = [0,ones(1,N1-1)];
c2 = [0,ones(1,N2-1)];
c3 = [0,ones(1,N3-1)];
for n1 = 1:N1
    for n2 = 1:N2
        for n3 = 1:N3
            for j3 = 1:N3
                n3minusj3 = wrapDFT(n3-j3,N3)+1+N3/2;
                part1 = -2*pi*1i*(n1-1-N1/2)*c1(n1)*v/L1*ycos(n3minusj3);
                part2 = -2*pi*1i*(n2-1-N2/2)*c2(n2)*v/L2*ysin(n3minusj3);
                A(n3,j3,n1,n2) = part1+part2;
            end
            
            part3 = -1i*(n3-1-N3/2)*c3(n3)*u(mode);
            part4 = -0.5*sigma^2*(n3-1-N3/2)^2;
            A(n3,n3,n1,n2) = A(n3,n3,n1,n2)+part3+part4;
        end
    end
end

% exponential of the coefficient matrix
expA = zeros(N3,N3,N1,N2);
for n1 = 1:N1
    for n2 = 1:N2
        expA(:,:,n1,n2) = expm(A(:,:,n1,n2)*Lt/(Nt-1));
    end
end
clear A;

% propagation
for nt = 2:Nt
    for n1 = 1:N1
        for n2 = 1:N2
            y(n1,n2,:,nt) = reshape(expA(:,:,n1,n2)*reshape(y(n1,n2,:,nt-1),[],1),1,1,[]);
        end
    end
end

% Inverse Fourier transform
for nt = 2:Nt
    fx(:,:,:,nt) = ifftn(ifftshift(y(:,:,:,nt)./shift1./shift2./shift3*n1*n2*n3),'symmetric');
end

% plot
for nt = 1:Nt
    figure;
    plot(x3,reshape(sum(sum(fx(:,:,:,nt)*L1/N1*L2/N2,1),2),[],1,1));
end

for nt = 1:Nt
    figure;
    surf(x1,x2,sum(fx(:,:,:,nt)*2*pi/N3,3)');
    view([0,0,1]);
end

rmpath('..\lib');

end

