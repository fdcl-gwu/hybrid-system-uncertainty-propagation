function [ xEst, fx, tTot, tIte ] = estimateMC( xMea )
close all;
addpath('..','..\..\lib');

timerTot = tic;

p = getParameter(1);
% parameters
v = p.v;
u = p.u;
sigma = p.sigma;
xo1 = p.xo1;
xo2 = p.xo2;
No = length(xo1);
xL1 = p.xL1;
xL2 = p.xL2;
sigmaL = p.sigmaL;
kL = p.kL;
nSample = p.nSample;

% grid
N1 = p.N1; N2 = p.N2;
L1 = p.L1; L2 = p.L2;
x1 = linspace(-L1/2,L1/2-L1/N1,N1);
x2 = linspace(-L2/2,L2/2-L2/N2,N2);
N3 = p.N3;
x3 = linspace(-pi,pi-2*pi/N3,N3);
Nt = p.Nt;
Lt = p.Lt;
t = linspace(0,Lt,Nt); dt = Lt/(Nt-1);

% likelihood function and true measurement
l = @(d,alpha,x1,x2) (1/sqrt(2*pi)/sigmaL)*exp(-(sqrt((x1-xL1).^2+(x2-xL2).^2)-d).^2/2/sigmaL^2) * ...
    (1/(2*pi*besseli(0,kL))).*exp(kL*cos(alpha-atan2(x2-xL2,x1-xL1)));

% draw initial samples
x = zeros(nSample,4,Nt);
x(:,:,1) = [rand(nSample,1)*L1-L1/2,rand(nSample,1)*L1-L1/2,rand(nSample,1)*2*pi,randi(3,nSample,1)];

% pre-allocate memory
fx = zeros(N1,N2,N3,3,Nt);
fx(:,:,:,:,1) = 1/L1/L2/(2*pi)/3;
xEst = zeros(Nt,4);
xEst(1,4) = 1;
tIte = zeros(Nt-1,1);

% propagation and estimation
for nt = 2:Nt
    timerIte = tic;
    
    % continuous propagation
    Bt = randn(nSample,1);
    for ns = 1:nSample
        x(ns,3,nt) = x(ns,3,nt-1) + u(x(ns,4,nt-1))*dt + Bt(ns)*sigma*sqrt(dt);
    end
    x(:,1,nt) = x(:,1,nt-1) + v*dt*(sin(x(:,3,nt))-sin(x(:,3,nt-1)))./(x(:,3,nt)-x(:,3,nt-1));
    x(:,2,nt) = x(:,2,nt-1) - v*dt*(cos(x(:,3,nt))-cos(x(:,3,nt-1)))./(x(:,3,nt)-x(:,3,nt-1));
    
    % discrete propagation
    theta = atan2(xo2-x(:,2,nt),xo1-x(:,1,nt));
    dtheta = wrapToPi(theta-x(:,3,nt));

    lamda = getLamda(x(:,1,nt),x(:,2,nt),xo1,xo2);
    [ind2,~] = ind2sub([4,nSample],find(lamda'));
    uni = rand(nSample,1);
    
    x(:,4,nt) = x(:,4,nt-1);
    
    indInPotent = find((x(:,4,nt-1)==1 & ind2~=4));
    indIn = uni(indInPotent)>exp(-lamda(sub2ind([nSample,4],indInPotent,ind2(indInPotent)))*dt);
    indIn = indInPotent(indIn);
    ind1To2 = dtheta(sub2ind([nSample,3],indIn,ind2(indIn)))<0 & dtheta(sub2ind([nSample,3],indIn,ind2(indIn)))<0;
    x(indIn(ind1To2),4,nt) = 2;
    x(indIn(~ind1To2),4,nt) = 3;
    
    indOutPotent = find((x(:,4,nt-1)~=1 & ind2==4));
    indOut = uni(indOutPotent)>exp(-lamda(indOutPotent,4)*dt);
    indOut = indOutPotent(indOut);
    x(indOut,4,nt) = 1;
    
    % measurement update
    w = ones(nSample,1)/nSample.*l(xMea(nt,1),xMea(nt,2),x(:,1,nt),x(:,2,nt));
    w = w/sum(w);
    
    % density
    for ns = 1:nSample
        [~,index1] = min(abs(x(ns,1,nt)-x1));
        [~,index2] = min(abs(x(ns,2,nt)-x2));
        [~,index3] = min(abs(wrapToPi(x(ns,3,nt)-x3)));
        fx(index1,index2,index3,x(ns,4,nt),nt) = fx(index1,index2,index3,x(ns,4,nt),nt)+w(ns);
    end
    fx(:,:,:,:,nt) = fx(:,:,:,:,nt)/(L1/N1*L2/N2*(2*pi)/N3);
    
    % estimation
    xEst(nt,1) = sum(x(:,1,nt).*w);
    xEst(nt,2) = sum(x(:,2,nt).*w);
    xEst(nt,3) = atan2(sum(sin(x(:,3,nt)).*w),sum(cos(x(:,3,nt)).*w));
    [~,xEst(nt,4)] = max(sum([x(:,4,nt)==1,x(:,4,nt)==2,x(:,4,nt)==3]));
    
    % re-sampling
    x(:,:,nt) = resample(x(:,:,nt),w,nSample);
    
    tIte(nt-1) = toc(timerIte);
end

tTot = toc(timerTot);

rmpath('..','..\..\lib');

end

