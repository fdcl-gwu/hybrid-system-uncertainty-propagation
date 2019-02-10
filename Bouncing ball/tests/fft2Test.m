clear;
close all;

% Gaussian distribution parameters
miu = [0;0];
sigma = [1,0.5;0.5,1];

% density function
f = @(x)1/(2*pi)/sqrt(det(sigma))*exp(-1/2*(x-miu)'*sigma^-1*(x-miu));

% sampling
n1 = 100; n2 = 100;
N1 = ceil(n1/2)-1; N2 = ceil(n2/2)-1;
L1 = 10; L2 = 10;
x1 = linspace(-L1/2,L1/2,n1);
x2 = linspace(-L2/2,L2/2,n2);
fx = zeros(n1,n2);
for i = 1:n1
    for j = 1:n2
        fx(i,j) = f([x1(i);x2(j)]);
    end
end

% 2D fft
y = fftshift(fft2(fx))/n1/n2;

% reconstruct density function
if mod(n1,2) == 0
    y(1,:) = [];
end
fraq1 = (-N1:N1)'/L1*2*N1/(2*N1+1);

if mod(n2,2) == 0
    y(:,1) = [];
end
fraq2 = (-N2:N2)/L2*2*N2/(2*N2+1);

f_fft = @(x)sum(sum(y.*exp(1i*fraq1*2*pi*x(1)).*exp(1i*fraq2*2*pi*x(2))));
f_fftx = zeros(n1,n2);
for i = 1:n1
    for j = 1:n2
        f_fftx(i,j) = real(f_fft([x1(i)+L1/2;x2(j)+L2/2]));
    end
end

% plot
figure; hold on;
surf(x1,x2,fx,'FaceColor','b','FaceAlpha',0.5);
surf(x1,x2,f_fftx,'FaceColor','r','FaceAlpha',0.5);

