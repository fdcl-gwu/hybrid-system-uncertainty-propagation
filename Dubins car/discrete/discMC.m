function [ p, x ] = discMC(  )
close all;
rng('shuffle');
addpath('..\..\lib');

% parameters
xo1 = [0,-1.5,1.5];
xo2 = [0.5,1,1];
No = length(xo1);
d = 0.5;
nSample = 100000;

% grid
N1 = 100; N2 = 100;
L1 = 6; L2 = 6;
x1 = linspace(-L1/2,L1/2-L1/N1,N1);
x2 = linspace(-L2/2,L2/2-L2/N2,N2);
N3 = 50;
x3 = linspace(-pi,pi-2*pi/N3,N3);

% initial conditions
x1_0 = 0; x2_0 = 0;
sigma1_0 = 0.2; sigma2_0 = 0.2;
x3_0 = pi/2;
k_0 = 20;
s_0 = 1;

% draw samples from initial condition
x = zeros(nSample,4,2);
x(:,1,1) = normrnd(x1_0,sigma1_0,nSample,1);
x(:,2,1) = normrnd(x2_0,sigma2_0,nSample,1);
x(:,3,1) = vmrnd(x3_0,k_0,nSample);
x(:,4,1) = ones(nSample,1)*s_0;

% propagate samples
theta = atan2(xo2-x(:,2,1),xo1-x(:,1,1));
dtheta = wrapToPi(theta-x(:,3,1));
in = false(nSample,No);
for no = 1:No
    in(:,no) = sqrt(sum((x(:,1:2,1)-[xo1(no),xo2(no)]).^2,2)) < d;
end

mode1 = x(:,4,1) == 1;
mode2 = x(:,4,1) == 2;
mode3 = x(:,4,1) == 3;

x(:,4,2) = x(:,4,1);
x(~sum(in,2) & (mode2 | mode3),4,2) = 1;

[Ind1,~] = find(in & mode1);
x(Ind1(dtheta(in & mode1)<0 & dtheta(in & mode1)>=-pi),4,2) = 2;
x(Ind1(dtheta(in & mode1)>=0 & dtheta(in & mode1)<pi),4,2) = 3;

% probability of each mode
p(1) = nnz(x(:,4,2)==1)/nSample;
p(2) = nnz(x(:,4,2)==2)/nSample;
p(3) = nnz(x(:,4,2)==3)/nSample;

rmpath('..\..\lib');

end

