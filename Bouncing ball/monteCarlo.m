clear; close all;

% parameters
g = 9.8;                                    % Gravity constant
niu = 0.05;                                 % Air drag coefficient
sigmaNiu = 0.01;                            % standard deviation of air drag coefficient
c = 0.9;                                    % coefficient of restitution
sigmaV = 0.5;                               % standard deviation for velocity reset
x0 = [1.5;0];                               % initial condition
sigma0 = [0.2^2,0;0,0.5^2];                 % covariance matrix of initial condition
nSample = 1000000;                           % sample size

% grid
n1 = 100; n2 = 100;
L1 = 5; L2 = 16;
x1 = linspace(-L1/2,L1/2-L1/n1,n1).'; x1(abs(x1)<1e-10) = 0;
x2 = linspace(-L2/2,L2/2-L2/n2,n2); x2(abs(x2)<1e-10) = 0;
nt = 241;
Lt = 6; dt = Lt/(nt-1);
t = linspace(0,Lt,nt);

% state variable
x = zeros(nSample,2,nt);

% draw samples from initial density
x(:,:,1) = mvnrnd(x0,sigma0,nSample);

% sample propagation
rng('shuffle');
for i = 2:nt
    x(:,1,i) = x(:,1,i-1)+x(:,2,i-1)*dt;
    
    niu_s = normrnd(niu,sigmaNiu,nSample,1);
    x(:,2,i) = x(:,2,i-1)-g*dt-niu_s.*x(:,2,i-1).*abs(x(:,2,i-1))*dt;
    
    index = find(x(:,1,i)<0);
    x(index,1,i) = -x(index,1,i);
    x(index,2,i) = normrnd(-c*x(index,2,i),sigmaV*ones(length(index),1),length(index),1);
end

% density approximation
fx = zeros(n1,n2,nt);
for i = 1:nt
    fx1 = x(:,1,i)>repmat(x1',nSample,1) & x(:,1,i)<=repmat(x1'+L1/n1,nSample,1);
    fx2 = x(:,2,i)>repmat(x2,nSample,1) & x(:,2,i)<=repmat(x2+L2/n2,nSample,1);
    parfor j = 1:n1
        fx(j,:,i) = sum(fx1(:,j) & fx2)/nSample*n1*n2/L1/L2;
    end
    fprintf(strcat(num2str(i),'th iteration finished\n'));
end

% plot
for i = 1:nt
    figure;
    surf(x2,x1,fx(:,:,i));
    view([0,0,1]);
end


