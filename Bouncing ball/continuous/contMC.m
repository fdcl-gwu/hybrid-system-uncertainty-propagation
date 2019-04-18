function [fx] = contMC()

close all;
rng('shuffle');
addpath('..');

p = getParameter(1);
% parameters
g = p.g;
niu = p.niu;
sigma = p.sigma;
x0 = p.x0;
sigma0 = p.sigma0;
nSample = p.nSample;

% grid
N1 = p.N1; N2 = p.N2;
L1 = p.L1; L2 = p.L2;
x1 = linspace(-L1/2,L1/2-L1/N1,N1)';
x2 = linspace(-L2/2,L2/2-L2/N2,N2);
Nt = 41;
Lt = 1; dt = Lt/(Nt-1);
t = linspace(0,Lt,Nt);

% sample from initial density
x = zeros(nSample,2,Nt);
x(:,:,1) = mvnrnd(x0,sigma0,nSample);

% sample propagation
for nt = 2:Nt
    Bt = randn(nSample,1);
    x(:,2,nt) = x(:,2,nt-1) - (g+niu*x(:,2,nt-1).*abs(x(:,2,nt-1)))*dt + sigma*x(:,2,nt-1).^2.*Bt*sqrt(dt);
    x(:,1,nt) = x(:,1,nt-1) + (x(:,2,nt)+x(:,2,nt-1))/2*dt;
end

% density approximation
fx = zeros(N1,N2,nt);
for nt = 1:Nt
    for ns = 1:nSample
        [~,index1] = min(abs(x(ns,1,nt)-x1));
        [~,index2] = min(abs(x(ns,2,nt)-x2));

        fx(index1,index2,nt) = fx(index1,index2,nt)+1;
    end
    fx(:,:,nt) = fx(:,:,nt)/nSample*N1*N2/L1/L2;
    fprintf(strcat(num2str(nt),'th iteration finished\n'));
end

% plot
for nt = 1:Nt
    figure;
    surf(x2,x1,fx(:,:,nt));
    view([0,0,1]);
end

rmpath('..');

end

