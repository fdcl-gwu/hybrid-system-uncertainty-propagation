function [fx] = contMC()

close all;
rng('shuffle');

% parameters
g = 9.8;
niu = 0.05;
sigma = 0.01;
x0 = [1.5;0];
sigma0 = [0.2^2,0;0,0.5^2];
nSample = 100000;

% grid
n1 = 100; n2 = 100;
L1 = 5; L2 = 16;
x1 = linspace(-L1/2,L1/2-L1/n1,n1)';
x2 = linspace(-L2/2,L2/2-L2/n2,n2);
nt = 41;
Lt = 1; dt = Lt/(nt-1);
t = linspace(0,Lt,nt);

% sample from initial density
x = zeros(nSample,2,nt);
x(:,:,1) = mvnrnd(x0,sigma0,nSample);

% sample propagation
for i = 2:nt
    Bt = randn(nSample,1);
    x(:,2,i) = x(:,2,i-1) - (g+niu*x(:,2,i-1).*abs(x(:,2,i-1)))*dt + sigma*x(:,2,i-1).^2.*Bt*sqrt(dt);
    x(:,1,i) = x(:,1,i-1) + (x(:,2,i)+x(:,2,i-1))/2*dt;
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

% plot
for i = 1:nt
figure;
surf(x2,x1,fx(:,:,i));
view([0,0,1]);
end

end

