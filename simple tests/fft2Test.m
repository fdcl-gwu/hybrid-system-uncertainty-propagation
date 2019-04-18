clear;
close all;

% Gaussian distribution parameters
miu = [2;2];
sigma = [1,0.5;0.5,1];

% density function
f = @(x)1/(2*pi)/sqrt(det(sigma))*exp(-1/2*(x-miu)'*sigma^-1*(x-miu));

% sampling
n1 = 200; n2 = 100;
L1 = 10; L2 = 10;
x1 = linspace(-L1/2,L1/2-L1/n1,n1);
x2 = linspace(-L2/2,L2/2-L2/n2,n2);
fx = zeros(n1,n2);
for i = 1:n1
    for j = 1:n2
        fx(i,j) = f([x1(i);x2(j)]);
    end
end

% 2D fft
shift1 = (-1).^((0:n1-1)-floor(n1/2)).';
shift2 = (-1).^((0:n2-1)-floor(n2/2));
y = fftshift(fft2(fx))/n1/n2;
y = shift1.*shift2.*y;

% reconstruct density function
freq1 = ((0:n1-1)'-floor(n1/2))/L1;
freq2 = ((0:n2-1)-floor(n2/2))/L2;

f_fft = @(x)sum(sum(y.*exp(1i*freq1*2*pi*x(1)).*exp(1i*freq2*2*pi*x(2))));
f_fftx = zeros(n1,n2);
for i = 1:n1
    for j = 1:n2
        f_fftx(i,j) = real(f_fft([x1(i);x2(j)]));
    end
end

% plot
figure; hold on;
surf(x2,x1,fx,'FaceColor','b','FaceAlpha',0.5);
surf(x2,x1,f_fftx,'FaceColor','r','FaceAlpha',0.5);

