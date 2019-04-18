function [ fx ] = discreteMC(  )

close all;
addpath('..');

p = getParameter(1);
% parameters
c = p.c;                                    % coefficient of restitution
sigmaV = p.sigmaV;                          % standard deviation of coefficient of restitution
x0 = [0;-3];                                % initial condition
sigma0 = [0.3^2,0.04;0.04,0.5^2];           % covariance matrix of initial condition
nSample = 100000;                           % sample size

% grid
N1= p.N1; N2 = p.N2;
L1 = p.L1; L2 = p.L2;
x1 = linspace(-L1/2,L1/2-L1/N1,N1).'; x1(abs(x1)<1e-10) = 0;
x2 = linspace(-L2/2,L2/2-L2/N2,N2); x2(abs(x2)<1e-10) = 0;
dt = 1/40;

% draw samples from initial density
x(:,:,1) = mvnrnd(x0,sigma0,nSample);

% sample propagation
index = x(:,1,1)<0;
x(index,1,2) = -x(index,1,1);
x(index,2,2) = normrnd(-c*x(index,2,1),sigmaV*ones(sum(index),1));
x(~index,:,2) = x(~index,:,1);

% density approximation
fx = zeros(N1,N2,2);
for nt = 1:2
    for ns = 1:nSample
        [~,index1] = min(abs(x(ns,1,nt)-x1));
        [~,index2] = min(abs(x(ns,2,nt)-x2));
        
        fx(index1,index2,nt) = fx(index1,index2,nt)+1;
    end
    fx(:,:,nt) = fx(:,:,nt)/nSample*N1*N2/L1/L2;
end

% plot
figure;
surf(x2,x1,fx(:,:,1));
view([0,0,1]);

figure;
surf(x2,x1,fx(:,:,2));
view([0,0,1]);

rmpath('..');

end

