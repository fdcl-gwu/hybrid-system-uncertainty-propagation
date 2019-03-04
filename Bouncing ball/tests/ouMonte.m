clear; close all;
rng('shuffle');

% parameters
miu = -2;
theta = 0.5;
sigma = 1;
x0 = 1;
D = sigma^2/2;
nSample = 1000000;

% density function
f = @(t,x)sqrt(theta/(2*pi*D*(1-exp(-2*theta*t))))*exp(-theta/(2*D)*...
    (x-miu-(x0-miu)*exp(-theta*t)).^2/(1-exp(-2*theta*t)));

% grid
Lx = 10;
nx = 100;
x = linspace(-Lx/2,Lx/2-Lx/nx,nx)';
Lt = 1;
nt = 40; dt = Lt/nt;
t = linspace(Lt/nt,Lt,nt);

% true density values
fxTure = zeros(nx,nt);
for i = 1:nt
    fxTure(:,i) = f(t(i),x);
end

% sampling from initial density
s = zeros(nSample,nt);
s(:,1) = randn(nSample,1)*sqrt(D/theta*(1-exp(-2*theta*t(1))))+miu+(x0-miu)*exp(-theta*t(1));

% sample propagation
for i = 2:nt
    Wt = randn(nSample,1);
    s(:,i) = s(:,i-1) + theta*(miu-s(:,i-1))*dt + sigma*Wt*sqrt(dt);
end

% reconstruct density
fx = zeros(nx,nt);
for i = 1:nt
    for j = 1:nx
        fx(j,i) = sum(s(:,i)>=x(j)-Lx/nx/2 & s(:,i)<x(j)+Lx/nx/2)/nSample*nx/Lx;
    end
end

% plot
figure; hold on;
for i = 1:nt
    plot3(x,ones(nx,1)*t(i),fxTure(:,i),'b');
    plot3(x,ones(nx,1)*t(i),fx(:,i),'r');
end

