function [ fx ] = discreteMC(  )

close all;

% parameters
c = 0.95;                                   % coefficient of restitution
sigmaV = 0.5;                               % standard deviation of coefficient of restitution
x0 = [0;-3];                                % initial condition
sigma0 = [0.3^2,0.04;0.04,0.5^2];           % covariance matrix of initial condition
nSample = 100000;                           % sample size

% grid
n1 = 100; n2 = 50;
L1 = 5; L2 = 16;
x1 = linspace(-L1/2,L1/2-L1/n1,n1).'; x1(abs(x1)<1e-10) = 0;
x2 = linspace(-L2/2,L2/2-L2/n2,n2); x2(abs(x2)<1e-10) = 0;
dt = 1/40;

% draw samples from initial density
x(:,:,1) = mvnrnd(x0,sigma0,nSample);

% sample propagation
index = x(:,1,1)<0;
x(index,1,2) = -x(index,1,1);
x(index,2,2) = normrnd(-c*x(index,2,1),sigmaV*ones(sum(index),1));
x(~index,:,2) = x(~index,:,1);

% density approximation
fx = zeros(n1,n2,2);
for i = 1:2
    fx1 = x(:,1,i)>repmat(x1'-L1/n1/2,nSample,1) & x(:,1,i)<=repmat(x1'+L1/n1/2,nSample,1);
    fx2 = x(:,2,i)>repmat(x2-L2/n2/2,nSample,1) & x(:,2,i)<=repmat(x2+L2/n2/2,nSample,1);
    parfor j = 1:n1
        fx(j,:,i) = sum(fx1(:,j) & fx2)/nSample*n1*n2/L1/L2;
    end
end

% plot
figure;
surf(x2,x1,fx(:,:,1));
view([0,0,1]);

figure;
surf(x2,x1,fx(:,:,2));
view([0,0,1]);

end

