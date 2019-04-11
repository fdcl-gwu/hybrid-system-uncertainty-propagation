function [ fx, x ] = hybridMC(  )

close all;
rng('shuffle');
addpath('..\..\lib');

% parameters
v = 1;
u = [0,2,-2];
sigma = 0.2;
xo1 = [0,-1.5,1.5];
xo2 = [0,1,1];
No = length(xo1);
epsilonIn = 6;
cIn = 400;
epsilonOut = 0.3;
cOut = 100;
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
x = zeros(nSample,4,Nt);
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

    distance = zeros(nSample,No);
    lamdaIn = zeros(nSample,3);
    lamdaOut = cOut*ones(nSample,1);
    for no = 1:No
        distance(:,no) = sqrt(sum((x(:,1:2,nt)-[xo1(no),xo2(no)]).^2,2));
        lamdaIn(:,no) = cIn*exp(-distance(:,no)*epsilonIn);
        lamdaOut = lamdaOut.*exp(-epsilonOut./distance(:,no));
    end

    % propagate samples
    in(:,1) = poissrnd(lamdaIn(:,1)*dt);
    in(:,2) = poissrnd(lamdaIn(:,2)*dt);
    in(:,3) = poissrnd(lamdaIn(:,3)*dt);
    out = poissrnd(lamdaOut*dt);

    x(:,4,nt) = x(:,4,nt-1);
    parfor ns = 1:nSample
        timeIn1 = rand(1,in(ns,1));
        timeIn2 = rand(1,in(ns,2));
        timeIn3 = rand(1,in(ns,3));
        timeOut = rand(1,out(ns));
        time = [timeIn1,timeIn2,timeIn3,timeOut
            1*ones(1,length(timeIn1)),2*ones(1,length(timeIn2)),3*ones(1,length(timeIn3)),4*ones(1,length(timeOut))];
        [~,ind] = sort(time(1,:));
        time = time(:,ind);

        if isempty(timeOut)
            if x(ns,4,nt) == 1 && ~isempty(time)
                if dtheta(ns,time(2,1))<0 && dtheta(ns,time(2,1))>=-pi
                    x(ns,4,nt) = 2;
                else
                    x(ns,4,nt) = 3;
                end
            end
        else
            if time(2,end) == 4
                x(ns,4,nt) = 1;
            else
                lastOutIndex = find(time(2,:)==4,1,'last');
                if dtheta(ns,time(2,lastOutIndex+1))<0 && dtheta(ns,time(2,lastOutIndex+1))>=-pi
                    x(ns,4,nt) = 2;
                else
                    x(ns,4,nt) = 3;
                end
            end
        end
    end
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
        plot3(xo1(no)+0.5*cos(0:0.01:2*pi),xo2(no)+0.5*sin(0:0.01:2*pi),ones(1,length(0:0.01:2*pi)),'Color','k','LineWidth',3);
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
save(strcat('D:\result-dubins car\',sprintf('%i-%i-%i-%i-%i-%i',round(clock)),'-MC','.mat'),'parameter','fx');

rmpath('..\..\lib');

end

