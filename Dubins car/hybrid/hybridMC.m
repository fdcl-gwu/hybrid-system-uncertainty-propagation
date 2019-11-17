function [ fx, x ] = hybridMC(  )

close all;
rng(1);
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
t = linspace(0,Lt,Nt);
dt = Lt/(Nt-1);

% initial conditions
x1_0 = p.x1_0; x2_0 = p.x2_0;
sigma1_0 = p.sigma1_0; sigma2_0 = p.sigma2_0;
x3_0 = p.x3_0;
k_0 = p.k_0;
s_0 = p.s_0;

% draw samples from initial condition
x = zeros(nSample,4,Nt);
x(:,1,1) = normrnd(x1_0,sigma1_0,nSample,1);
x(:,2,1) = normrnd(x2_0,sigma2_0,nSample,1);
x(:,3,1) = vmrnd(x3_0,k_0,nSample);
x(:,4,1) = ones(nSample,1)*s_0;

% pre-allocate memory
fx = zeros(N1,N2,N3,3,Nt);
tIte = zeros(Nt-1,1);

% recover initial density
for ns = 1:nSample
    [~,index1] = min(abs(x(ns,1,1)-x1));
    [~,index2] = min(abs(x(ns,2,1)-x2));
    [~,index3] = min(abs(wrapToPi(x(ns,3,1)-x3)));

    fx(index1,index2,index3,x(ns,4,1),1) = fx(index1,index2,index3,x(ns,4,1),1)+1;
end
fx(:,:,:,:,1) = fx(:,:,:,:,1)/nSample*N1/L1*N2/L2*N3/(2*pi);

% propagate samples
for nt = 2:Nt
    timerIte = tic;
    
    % continuous
    Bt = randn(nSample,1);
    for ns = 1:nSample
        x(ns,3,nt) = x(ns,3,nt-1) + u(x(ns,4,nt-1))*dt + Bt(ns)*sigma*sqrt(dt);
    end
    x(:,1,nt) = x(:,1,nt-1) + v*dt*(sin(x(:,3,nt))-sin(x(:,3,nt-1)))./(x(:,3,nt)-x(:,3,nt-1));
    x(:,2,nt) = x(:,2,nt-1) - v*dt*(cos(x(:,3,nt))-cos(x(:,3,nt-1)))./(x(:,3,nt)-x(:,3,nt-1));
    
    % discrete
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
    
    % recover density
    for ns = 1:nSample
        [~,index1] = min(abs(x(ns,1,nt)-x1));
        [~,index2] = min(abs(x(ns,2,nt)-x2));
        [~,index3] = min(abs(wrapToPi(x(ns,3,nt)-x3)));
        
        fx(index1,index2,index3,x(ns,4,nt),nt) = fx(index1,index2,index3,x(ns,4,nt),nt)+1;
    end
    fx(:,:,:,:,nt) = fx(:,:,:,:,nt)/nSample*N1/L1*N2/L2*N3/(2*pi);
    
    tIte(nt-1) = toc(timerIte);
end

tTot = toc(timerTot);

% plot
for nt = 1:4:Nt
    figure;
    plot(x3,reshape(sum(sum(sum(fx(:,:,:,:,nt)*L1/N1*L2/N2,1),2),4),[],1,1));
end

for nt = 1:4:Nt
    figure; hold on;
    for no = 1:No
        scatter3(xo1(no),xo2(no),1,'Marker','o','SizeData',20,'MarkerFaceColor','k','MarkerEdgeColor','k');
        plot3(xo1(no)+0.5*cos(0:0.01:2*pi),xo2(no)+0.5*sin(0:0.01:2*pi),ones(1,length(0:0.01:2*pi)),'Color','k','LineWidth',3);
    end
    surf(x1,x2,sum(sum(fx(:,:,:,:,nt)*2*pi/N3,3),4)');
    view([0,0,1]);
end

% save data
save(strcat('D:\result-dubins car\',sprintf('%i-%i-%i-%i-%i-%i',round(clock)),'-MC','.mat'),'p','fx','tTot','tIte');

rmpath('..','..\..\lib');

end

