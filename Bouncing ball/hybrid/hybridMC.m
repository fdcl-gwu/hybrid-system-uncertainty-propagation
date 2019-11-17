function [fx] = hybridMC()

close all;
rng(1);
addpath('..');
timerTot = tic;

p = getParameter(1);
% parameters
g = p.g;                                    % Gravity constant
niu = p.niu;                                % Air drag coefficient
sigmaNiu = p.sigma;                         % standard deviation of air drag coefficient
c = p.c;                                    % coefficient of restitution
sigmaV = p.sigmaV;                          % standard deviation for velocity reset
x0 = p.x0;                                  % initial condition
sigma0 = p.sigma0;                          % covariance matrix of initial condition
nSample = p.nSample;                        % sample size

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

% pre-allocate memory
fx = zeros(N1,N2,Nt);
tIte = zeros(Nt-1,1);

% recover initial density
for ns = 1:nSample
    [~,index1] = min(abs(x(ns,1,1)-x1));
    [~,index2] = min(abs(x(ns,2,1)-x2));

    fx(index1,index2,1) = fx(index1,index2,1)+1;
end
fx(:,:,1) = fx(:,:,1)/nSample*N1*N2/L1/L2;

% sample propagation
for nt = 2:Nt
    timerIte = tic;
    
    % continuous propagation
    Bt = randn(nSample,1);
    x(:,2,nt) = x(:,2,nt-1) - (g+niu*x(:,2,nt-1).*abs(x(:,2,nt-1)))*dt + sigmaNiu*x(:,2,nt-1).*abs(x(:,2,nt-1)).*Bt*sqrt(dt);
    x(:,1,nt) = x(:,1,nt-1) + (x(:,2,nt)+x(:,2,nt-1))/2*dt;
    
    % discrete propagation
    index = find(x(:,1,nt)<0);
    x(index,1,nt) = -x(index,1,nt);
    x(index,2,nt) = normrnd(-c*x(index,2,nt),sigmaV*ones(length(index),1));
    
    % recover density
    for ns = 1:nSample
        [~,index1] = min(abs(x(ns,1,nt)-x1));
        [~,index2] = min(abs(x(ns,2,nt)-x2));
        
        fx(index1,index2,nt) = fx(index1,index2,nt)+1;
    end
    fx(:,:,nt) = fx(:,:,nt)/nSample*N1*N2/L1/L2;
    
    tIte(nt-1) = toc(timerIte);
end

tTot = toc(timerTot);

% plot
for nt = 1:Nt
    figure;
    surf(x2,x1,fx(:,:,nt));
    view([0,0,1]);
end

% save data
save(strcat('D:\result-bouncing ball\',sprintf('%i-%i-%i-%i-%i-%i',round(clock)),'.mat'),'p','x','fx','tTot','tIte','-v7.3');

rmpath('..');

end
