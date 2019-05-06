function [ fx, xTrue, xEst ] = estimateMC(  )
close all;
rng(4);
addpath('..','..\..\lib');
tic;

p = getParameter(1);
% parameters
g = p.g;                                    % Gravity constant
niu = p.niu;                                % Air drag coefficient
sigmaNiu = p.sigma;                         % standard deviation of air drag coefficient
c = p.c;                                    % coefficient of restitution
sigmaV = p.sigmaV;                          % standard deviation for velocity reset
sigmaM = p.sigmaM;                          % standard deviation of measurement
nSample = p.nSample;                        % number of samples

% grid
N1 = p.N1; N2 = p.N2;
L1 = p.L1; L2 = p.L2;
x1 = linspace(-L1/2,L1/2-L1/N1,N1).'; x1(abs(x1)<1e-10) = 0;
x2 = linspace(-L2/2,L2/2-L2/N2,N2); x2(abs(x2)<1e-10) = 0;
Nt = p.Nt;
Lt = p.Lt;
t = linspace(0,Lt,Nt); dt = Lt/(Nt+1);

% true state
xTrue = generateSample(1,p);
xTrue = reshape(xTrue,2,Nt)';

% measurement
xMea = randn(Nt,1)*sigmaM+xTrue;

% likelihood function
l = @(h,x1) 1/sqrt(2*pi)/sigmaM*exp(-(x1-h).^2/2/sigmaM^2);

% draw samples from initial distribution
x = zeros(nSample,2,Nt);
x(:,:,1) = [rand(nSample,1)*L1/2,rand(nSample,1)*L2-L2/2];

% propagation and estimation
fx = zeros(N1,N2,Nt);
xEst = zeros(Nt,2);
for nt = 2:Nt
    % propagate samples
    Bt = randn(nSample,1);
    x(:,2,nt) = x(:,2,nt-1) - (g+niu*x(:,2,nt-1).*abs(x(:,2,nt-1)))*dt + sigmaNiu*x(:,2,nt-1).*abs(x(:,2,nt-1)).*Bt*sqrt(dt);
    x(:,1,nt) = x(:,1,nt-1) + (x(:,2,nt)+x(:,2,nt-1))/2*dt;
    
    index = find(x(:,1,nt)<0);
    x(index,1,nt) = -x(index,1,nt);
    x(index,2,nt) = normrnd(-c*x(index,2,nt),sigmaV*ones(length(index),1));
    
    % density
    for ns = 1:nSample
        [~,index1] = min(abs(x(ns,1,nt)-x1));
        [~,index2] = min(abs(x(ns,2,nt)-x2));
        fx(index1,index2,nt) = fx(index1,index2,nt)+1;
    end
    fx(:,:,nt) = fx(:,:,nt)/nSample*N1*N2/L1/L2;
    
    % measurement update
    lx = repmat(l(xMea(nt),x1),1,N2);
    fx(:,:,nt) = fx(:,:,nt).*lx;
    fx(:,:,nt) = fx(:,:,nt)/(sum(sum(fx(:,:,nt)*L1*L2/N1/N2)));
    
    % estimation
    [~,index1] = max(max(fx(:,:,nt),[],2),[],1);
    [~,index2] = max(max(fx(:,:,nt),[],1),[],2);
    xEst(nt,:) = [x1(index1),x2(index2)];
    
    % re-sampling
    x(:,:,nt) = randpdf2(x1,x2,fx(:,:,nt),nSample);
end

simulT = toc;

% plot
for nt = 1:4:Nt
    figure; hold on;
    surf(x2,x1,fx(:,:,nt));
    scatter3(xTrue(nt,2),xTrue(nt,1),max(max(fx(:,:,nt))),'MarkerFaceColor','r','MarkerEdgeColor','r');
    scatter3(xEst(nt,2),xEst(nt,1),max(max(fx(:,:,nt))),'MarkerFaceColor','b','MarkerEdgeColor','b');
    view([0,0,1]);
end

save(strcat('D:\result-bouncing ball\',sprintf('%i-%i-%i-%i-%i-%i',round(clock)),'-estimateMC.mat'),'fx','xTrue','xEst','xMea','p','simulT');

rmpath('..','..\..\lib');

end

