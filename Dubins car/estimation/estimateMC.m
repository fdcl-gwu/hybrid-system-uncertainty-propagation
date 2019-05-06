function [ fx, xTrue, xEst ] = estimateMC(  )
close all;
rng(10);
addpath('..','..\..\lib');
tic;

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

% true state
xTrue = generateSample(1);
xTrue = reshape(xTrue,4,Nt)';

% likelihood function and true measurement
l = @(d,alpha,x1,x2) (1/sqrt(2*pi)/sigmaL)*exp(-(sqrt((x1-xL1)^2+(x2-xL2)^2)-d)^2/2/sigmaL^2) * ...
    (1/(2*pi*besseli(0,kL)))*exp(kL*cos(alpha-atan2(x2-xL2,x1-xL1)));
xMea(:,1) = randn(Nt,1)*sigmaL+sqrt((xTrue(:,1)-xL1).^2+(xTrue(:,2)-xL2).^2);
for nt = 1:Nt
    xMea(nt,2) = vmrnd(atan2(xTrue(nt,2)-xL2,xTrue(nt,1)-xL1),kL,1);
end

% draw initial samples
x = zeros(nSample,4,Nt);
x(:,:,1) = [rand(nSample,1)*L1-L1/2,rand(nSample,1)*L1-L1/2,rand(nSample,1)*2*pi,randi(3,nSample,1)];

% propagation and estimation
fx = zeros(N1,N2,N3,3,Nt);
xEst = zeros(Nt,4);
for nt = 2:Nt
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
    lamdaIn = lamda(:,1:3);
    lamdaOut = lamda(:,4);

    in(:,1) = poissrnd(lamdaIn(:,1)*dt);
    in(:,2) = poissrnd(lamdaIn(:,2)*dt);
    in(:,3) = poissrnd(lamdaIn(:,3)*dt);
    out = poissrnd(lamdaOut*dt);

    x(:,4,nt) = x(:,4,nt-1);
    for ns = 1:nSample
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
    
    % density
    for ns = 1:nSample
        [~,index1] = min(abs(x(ns,1,nt)-x1));
        [~,index2] = min(abs(x(ns,2,nt)-x2));
        [~,index3] = min(abs(wrapToPi(x(ns,3,nt)-x3)));
        fx(index1,index2,index3,x(ns,4,nt),nt) = fx(index1,index2,index3,x(ns,4,nt),nt)+1;
    end
    fx(:,:,:,:,nt) = fx(:,:,:,:,nt)/nSample*N1/L1*N2/L2*N3/(2*pi);
    
    % measurement update
    lx = zeros(N1,N2,N3,3);
    for n1 = 1:N1
        for n2 = 1:N2
            lx(n1,n2,:,:) = l(xMea(nt,1),xMea(nt,2),x1(n1),x2(n2));
        end
    end
    
    fx(:,:,:,:,nt) = fx(:,:,:,:,nt).*lx;
    fx(:,:,:,:,nt) = fx(:,:,:,:,nt)/(sum(sum(sum(sum(fx(:,:,:,:,nt)))))*L1/N1*L2/N2*(2*pi)/N3);
    
    % estimation
    xEst(nt,1) = sum(x1'.*sum(sum(sum(fx(:,:,:,:,nt),4),3),2)*L2/N2*(2*pi)/N3)*L1/N1;
    xEst(nt,2) = sum(x2'.*reshape(sum(sum(sum(fx(:,:,:,:,nt),4),3),1),N2,1)*L1/N1*(2*pi)/N3)*L2/N2;
    xEst(nt,3) = atan2(sum(sin(x3)'.*reshape(sum(sum(sum(fx(:,:,:,:,nt),4),2),1),N3,1)*L1/N1*L2/N2)*(2*pi)/N3,...
        sum(cos(x3)'.*reshape(sum(sum(sum(fx(:,:,:,:,nt),4),2),1),N3,1)*L1/N1*L2/N2)*(2*pi)/N3);
    [~,xEst(nt,4)] = max(reshape(sum(sum(sum(fx(:,:,:,:,nt),3),2),1)*L1/N1*L2/N2*(2*pi)/N3,3,1));
    
    % re-sampling
    x(:,:,nt) = randpdfCar(x1,x2,x3,fx(:,:,:,:,nt),nSample);
end

simulT = toc;

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

save(strcat('D:\result-dubins car\',sprintf('%i-%i-%i-%i-%i-%i',round(clock)),'-estimateMC','.mat'),'fx','xTrue','xEst','p','simulT');

end

