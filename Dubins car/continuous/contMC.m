function [ fx, x ] = contMC( mode )
close all;
addpath('..','..\..\lib');
rng('shuffle');

p = getParameter(1);
% parameters
v = p.v;
u = p.u;
sigma = p.sigma;
nSample = p.nSample;
if ~exist('mode','var') || isempty('mode')
    mode = 1;
end

% grid
N1 = p.N1; N2 = p.N2;
L1 = p.L1; L2 = p.L2;
x1 = linspace(-L1/2,L1/2-L1/N1,N1);
x2 = linspace(-L2/2,L2/2-L2/N2,N2);
N3 = p.N3;
x3 = linspace(-pi,pi-2*pi/N3,N3);
Nt = 41;
Lt = 1;
t = linspace(0,Lt,Nt);
dt = Lt/(Nt-1);

% initial conditions
x1_0 = p.x1_0; x2_0 = p.x2_0;
sigma1_0 = p.sigma1_0; sigma2_0 = p.sigma2_0;
x3_0 = p.x3_0;
k_0 = p.k_0;

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

rmpath('..','..\..\lib');

end

