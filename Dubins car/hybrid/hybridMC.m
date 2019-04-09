function [ fx, x ] = hybridMC(  )

close all;
rng('shuffle');
addpath('..\..\lib');

% parameters
v = 1;
u = [0,2,-2];
sigma = 0.2;
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
Nt = 161;
Lt = 4;
t = linspace(0,Lt,Nt);
dt = Lt/(Nt-1);

% initial conditions
x1_0 = 0; x2_0 = -2;
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
for nt = 2:Nt
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
    in = false(nSample,No);
    for no = 1:No
        in(:,no) = sqrt(sum((x(:,1:2,nt)-[xo1(no),xo2(no)]).^2,2)) < d;
    end

    mode1 = x(:,4,nt-1) == 1;
    mode2 = x(:,4,nt-1) == 2;
    mode3 = x(:,4,nt-1) == 3;

    x(:,4,nt) = x(:,4,nt-1);
    x(~sum(in,2) & (mode2 | mode3),4,nt) = 1;
    
    [Ind1,~] = find(in & mode1);
    x(Ind1(dtheta(in & mode1)<0 & dtheta(in & mode1)>=-pi),4,nt) = 2;
    x(Ind1(dtheta(in & mode1)>=0 & dtheta(in & mode1)<pi),4,nt) = 3;
end

% convert sample to density
fx = zeros(N1,N2,N3,3,Nt);
for nt = 1:Nt
    for ns = 1:nSample
        [~,index1] = min(abs(x(ns,1,nt)-x1));
        [~,index2] = min(abs(x(ns,2,nt)-x2));
        [~,index3] = min(abs(wrapToPi(x(ns,3,nt)-x3)));
        
        fx(index1,index2,index3,x(ns,4,nt),nt) = fx(index1,index2,index3,x(ns,4,nt),nt)+1;
    end
    fx(:,:,:,:,nt) = fx(:,:,:,:,nt)/nSample*N1/L1*N2/L2*N3/(2*pi);
end

% plot
for nt = 1:4:Nt
    figure;
    plot(x3,reshape(sum(sum(sum(fx(:,:,:,:,nt)*L1/N1*L2/N2,1),2),4),[],1,1));
end

for nt = 1:4:Nt
    figure; hold on;
    for no = 1:No
        scatter3(xo1(no),xo2(no),1,'Marker','o','SizeData',20,'MarkerFaceColor','k','MarkerEdgeColor','k');
        plot3(xo1(no)+d*cos(0:0.01:2*pi),xo2(no)+d*sin(0:0.01:2*pi),ones(1,length(0:0.01:2*pi)),'Color','k','LineWidth',3);
    end
    surf(x1,x2,sum(sum(fx(:,:,:,:,nt)*2*pi/N3,3),4)');
    view([0,0,1]);
end

% save data
parameter.x1 = x1;
parameter.x2 = x2;
parameter.x3 = x3;
parameter.t = t;
parameter.xo1 = xo1;
parameter.xo2 = xo2;
parameter.d = d;
save(strcat('D:\result-dubins car\',sprintf('%i-%i-%i-%i-%i-%i',round(clock)),'-MC','.mat'),'parameter','fx');

rmpath('..\..\lib');

end

