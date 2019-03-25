function [ fx, x ] = contMC( mode )
close all;
addpath('..\..\lib');
rng('shuffle');

% parameters
v = 0.5;
u = [0,0.5,-0.5];
sigma = 0.2;
nSample = 100000;
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
dt = Lt/(Nt-1);
t = linspace(0,Lt,Nt);

% initial conditions
x1_0 = 0; x2_0 = -2;
sigma1_0 = 0.2; sigma2_0 = 0.2;
x3_0 = pi/2;
k_0 = 20;

% draw samples from initial condition
x = zeros(nSample,3,Nt);
x(:,1,1) = normrnd(x1_0,sigma1_0,nSample,1);
x(:,2,1) = normrnd(x2_0,sigma2_0,nSample,1);
x(:,3,1) = vmrnd(x3_0,k_0,nSample);

% propagate samples
for nt = 2:Nt
    Bt = randn(nSample,1);
    x(:,3,nt) = x(:,3,nt-1) + u(mode)*dt + Bt*sigma*sqrt(dt);
    x(:,1,nt) = x(:,1,nt-1) + v*dt*(sin(x(:,3,nt))-sin(x(:,3,nt-1)))./(x(:,3,nt)-x(:,3,nt-1));
    x(:,2,nt) = x(:,2,nt-1) - v*dt*(cos(x(:,3,nt))-cos(x(:,3,nt-1)))./(x(:,3,nt)-x(:,3,nt-1));
end

% convert sample to density
fx = zeros(N1,N2,N3,Nt);
for nt = 1:Nt
    for ns = 1:nSample
        [~,index1] = min(abs(x(ns,1,nt)-x1));
        [~,index2] = min(abs(x(ns,2,nt)-x2));
        [~,index3] = min(abs(wrapToPi(x(ns,3,nt)-x3)));
        
        fx(index1,index2,index3,nt) = fx(index1,index2,index3,nt)+1;
    end
    fx(:,:,:,nt) = fx(:,:,:,nt)/nSample*N1/L1*N2/L2*N3/(2*pi);
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

addpath('..\..\lib');

end

