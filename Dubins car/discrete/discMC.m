function [ prob, x ] = discMC(  )
% close all;
rng('shuffle');
addpath('..','..\..\lib');

p = getParameter(1);
% parameters
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
dt = 0.025;

% initial conditions
x1_0 = 0; x2_0 = -0.5;
sigma1_0 = p.sigma1_0; sigma2_0 = p.sigma2_0;
x3_0 = p.x3_0;
k_0 = p.k_0;
s_0 = 1;

% draw samples from initial condition
x = zeros(nSample,4,2);
x(:,1,1) = normrnd(x1_0,sigma1_0,nSample,1);
x(:,2,1) = normrnd(x2_0,sigma2_0,nSample,1);
x(:,3,1) = vmrnd(x3_0,k_0,nSample);
x(:,4,1) = ones(nSample,1)*s_0;

% theta
theta = atan2(xo2-x(:,2,1),xo1-x(:,1,1));
dtheta = wrapToPi(theta-x(:,3,1));

% transition rate
lamda = getLamda(x(:,1),x(:,2),xo1,xo2);
lamdaIn = lamda(:,1:3);
lamdaOut = lamda(:,4);

% propagate samples
in(:,1) = poissrnd(lamdaIn(:,1)*dt);
in(:,2) = poissrnd(lamdaIn(:,2)*dt);
in(:,3) = poissrnd(lamdaIn(:,3)*dt);
out = poissrnd(lamdaOut*dt);

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
        if x(ns,4,1) == 1 && ~isempty(time)
            if dtheta(ns,time(2,1))<0 && dtheta(ns,time(2,1))>=-pi
                x(ns,4,2) = 2;
            else
                x(ns,4,2) = 3;
            end
        else
            x(ns,4,2) = x(ns,4,1);
        end
    else
        if time(2,end) == 4
            x(ns,4,2) = 1;
        else
            lastOutIndex = find(time(2,:)==4,1,'last');
            if dtheta(ns,time(2,lastOutIndex+1))<0 && dtheta(ns,time(2,lastOutIndex+1))>=-pi
                x(ns,4,2) = 2;
            else
                x(ns,4,2) = 3;
            end
        end
    end
end

% probability of each mode
prob(1) = nnz(x(:,4,2)==1)/nSample;
prob(2) = nnz(x(:,4,2)==2)/nSample;
prob(3) = nnz(x(:,4,2)==3)/nSample;

rmpath('..','..\..\lib');

end

