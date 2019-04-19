function [ x ] = generateSample( n, p )

close all;

if ~exist('p','var')
    addpath('..');
    p = getParameter(1);
    rmpath('..');
end

% parameters
g = p.g;                                    % Gravity constant
niu = p.niu;                                % Air drag coefficient
sigmaNiu = p.sigma;                         % standard deviation of air drag coefficient
c = p.c;                                    % coefficient of restitution
sigmaV = p.sigmaV;                          % standard deviation for velocity reset
x0 = p.x0;                                  % initial condition
sigma0 = p.sigma0;                          % covariance matrix of initial condition
nSample = n;                                % sample size

% grid
N1 = p.N1; N2 = p.N2;
L1 = p.L1; L2 = p.L2;
x1 = linspace(-L1/2,L1/2-L1/N1,N1).'; x1(abs(x1)<1e-10) = 0;
x2 = linspace(-L2/2,L2/2-L2/N2,N2); x2(abs(x2)<1e-10) = 0;
Nt = p.Nt;
Lt = p.Lt; dt = Lt/(Nt-1);
t = linspace(0,Lt,Nt);

% state variable
x = zeros(nSample,2,Nt);

% draw samples from initial density
x(:,:,1) = mvnrnd(x0,sigma0,nSample);

% sample propagation
for i = 2:Nt
    Bt = randn(nSample,1);
    x(:,2,i) = x(:,2,i-1) - (g+niu*x(:,2,i-1).*abs(x(:,2,i-1)))*dt + sigmaNiu*x(:,2,i-1).*abs(x(:,2,i-1)).*Bt*sqrt(dt);
    x(:,1,i) = x(:,1,i-1) + (x(:,2,i)+x(:,2,i-1))/2*dt;
    
    index = find(x(:,1,i)<0);
    x(index,1,i) = -x(index,1,i);
    x(index,2,i) = normrnd(-c*x(index,2,i),sigmaV*ones(length(index),1));
end

end

