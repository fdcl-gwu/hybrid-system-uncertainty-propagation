function [fx] = hybridMC()

close all;
rng('shuffle');
tic;

% parameters
g = 9.8;                                    % Gravity constant
niu = 0.05;                                 % Air drag coefficient
sigmaNiu = 0.01;                            % standard deviation of air drag coefficient
c = 0.95;                                   % coefficient of restitution
sigmaV = 0.5;                               % standard deviation for velocity reset
x0 = [1.5;0];                               % initial condition
sigma0 = [0.2^2,0;0,0.5^2];                 % covariance matrix of initial condition
nSample = 1000000;                          % sample size

% grid
n1 = 100; n2 = 50;
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
for i = 2:nt
    Bt = randn(nSample,1);
    x(:,2,i) = x(:,2,i-1) - (g+niu*x(:,2,i-1).*abs(x(:,2,i-1)))*dt + sigmaNiu*x(:,2,i-1).*abs(x(:,2,i-1)).*Bt*sqrt(dt);
    x(:,1,i) = x(:,1,i-1) + (x(:,2,i)+x(:,2,i-1))/2*dt;
    
    index = find(x(:,1,i)<0);
    x(index,1,i) = -x(index,1,i);
    x(index,2,i) = normrnd(-c*x(index,2,i),sigmaV*ones(length(index),1));
end

% density approximation
fx = zeros(n1,n2,nt);
for i = 1:nt
    for j = 1:nSample
        [~,index1] = min(abs(x(j,1,i)-x1));
        [~,index2] = min(abs(x(j,2,i)-x2));
        
        fx(index1,index2,i) = fx(index1,index2,i)+1;
    end
    fx(:,:,i) = fx(:,:,i)/nSample*n1*n2/L1/L2;
    fprintf(strcat(num2str(i),'th iteration finished\n'));
end

simulT = toc;

% plot
for i = 1:nt
    figure;
    surf(x2,x1,fx(:,:,i));
    view([0,0,1]);
end

% save data
parameter.g = g;
parameter.niu = niu;
parameter.sigmaNiu = sigmaNiu;
parameter.c = c;
parameter.sigmaV = sigmaV;
parameter.x0 = x0;
parameter.sigma0 = sigma0;
parameter.nSample = nSample;

save(strcat('D:\result-bouncing ball\',sprintf('%i-%i-%i-%i-%i-%i',round(clock)),'.mat'),'parameter','x1','x2','t','x','fx','simulT','-v7.3');

end
