function [ x, xMea ] = generateSample( n, p )

addpath('..','..\..\lib');
close all;

if ~exist('p','var')
    p = getParameter(1);
end
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
nSample = n;

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

% initial measurement
xMea = zeros(nSample,2,Nt);
xMea(:,1,1) = sqrt((x(:,1,1)-xL1).^2+(x(:,2,1)-xL2).^2)+randn(nSample,1)*sigmaL;
for ns = 1:nSample
    xMea(ns,2,1) = vmrnd(atan2(x(ns,2,1)-xL2,x(ns,1,1)-xL1),kL,1);
end

% propagate samples
for nt = 2:Nt
    % continuous
    Bt = randn(nSample,1);
    for ns = 1:nSample
        x(ns,3,nt) = x(ns,3,nt-1) + u(x(ns,4,nt-1))*dt + Bt(ns)*sigma*sqrt(dt);
    end
    x(:,1,nt) = x(:,1,nt-1) + v*dt*(sin(x(:,3,nt))-sin(x(:,3,nt-1)))./(x(:,3,nt)-x(:,3,nt-1));
    x(:,2,nt) = x(:,2,nt-1) - v*dt*(cos(x(:,3,nt))-cos(x(:,3,nt-1)))./(x(:,3,nt)-x(:,3,nt-1));
    
    % theta
    theta = atan2(xo2-x(:,2,nt),xo1-x(:,1,nt));
    dtheta = wrapToPi(theta-x(:,3,nt));

    % rate function
    lamda = getLamda(x(:,1,nt),x(:,2,nt),xo1,xo2);
    lamdaIn = lamda(:,1:3);
    lamdaOut = lamda(:,4);

    % propagate samples
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
    
    % measurement
    xMea(:,1,nt) = sqrt((x(:,1,nt)-xL1).^2+(x(:,2,nt)-xL2).^2)+randn(nSample,1)*sigmaL;
    for ns = 1:nSample
        xMea(ns,2,nt) = vmrnd(atan2(x(ns,2,nt)-xL2,x(ns,1,nt)-xL1),kL,1);
    end
end

rmpath('..','..\..\lib');

end

