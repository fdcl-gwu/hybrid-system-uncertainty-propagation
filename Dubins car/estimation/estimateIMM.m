function [ xEst, miu, P, tTot, tIte ] = estimateIMM( xMea )
% close all;
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

% grid
N1 = p.N1; N2 = p.N2;
L1 = p.L1; L2 = p.L2;
x1 = linspace(-L1/2,L1/2-L1/N1,N1);
x2 = linspace(-L2/2,L2/2-L2/N2,N2);
N3 = p.N3;
x3 = linspace(-pi,pi-2*pi/N3,N3);
Nt = p.Nt;
Lt = p.Lt;
dt = Lt/(Nt-1);

% initial knowledge
miu = zeros(3,3,Nt);
miu(:,:,1) = repmat([0;0;0],1,3);
P = zeros(3,3,3,Nt);
P(:,:,:,1) = repmat([(L1*10)^2,0,0;0,(L2*10)^2,0;0,0,(20*pi)^2],1,1,3);
ps = zeros(3,Nt);
ps(:,1) = [1;1;1]/3;

% calculate lamda and kernel functions;
lamda = zeros(N1,N2,4);
for n2 = 1:N2
    lamda(:,n2,:) = reshape(getLamda(x1',ones(N2,1)*x2(n2),xo1,xo2),N1,1,4);
end
lamda = reshape(lamda,N1,N2,1,4);
kernel = getKernel(x1,x2,x3,xo1,xo2);

% Kalman filter parameters
Q = zeros(3,3); Q(3,3) = sigma^2*dt;
R = [sigmaL^2,0;0,getR(kL)];

% pre-allocate memory
xEst = zeros(4,Nt);
xEst(1,4) = 1;
tIte = zeros(Nt-1,1);

% IMM iteration
for nt = 2:Nt
    timerIte = tic;
    
    % propagate continuous dynamics
    for ns = 1:3
        miu(:,ns,nt) = miu(:,ns,nt-1)+[v*cos(miu(3,ns,nt-1));...
            v*sin(miu(3,ns,nt-1));u(ns)]*dt;
        F = [1,0,-v*sin(miu(3,ns,nt-1))*dt;
            0,1,v*cos(miu(3,ns,nt-1))*dt;
            0,0,1];
        P(:,:,ns,nt) = F*P(:,:,ns,nt-1)*F'+Q;
    end
    
    % propagate discrete dynamics
    PI = getPI(miu(:,:,nt),P(:,:,:,nt),lamda,kernel,dt,x1,x2,x3);
    ps(:,nt) = PI'*ps(:,nt-1);
    
    miuP = zeros(3,3);
    PP = zeros(3,3,3);
    for ns = 1:3
        miuP(:,ns) = (PI(1,ns)*ps(1,nt-1)*miu(:,1,nt)+...
            PI(2,ns)*ps(2,nt-1)*miu(:,2,nt)+...
            PI(3,ns)*ps(3,nt-1)*miu(:,3,nt))/ps(ns,nt);
        PP(:,:,ns) = (PI(1,ns)*ps(1,nt-1)*(P(:,:,1,nt)+(miu(:,1,nt)-miuP(:,ns))*(miu(:,1,nt)-miuP(:,ns))')+...
            PI(2,ns)*ps(2,nt-1)*(P(:,:,2,nt)+(miu(:,2,nt)-miuP(:,ns))*(miu(:,2,nt)-miuP(:,ns))')+...
            PI(3,ns)*ps(3,nt-1)*(P(:,:,3,nt)+(miu(:,3,nt)-miuP(:,ns))*(miu(:,3,nt)-miuP(:,ns))'))/ps(ns,nt);
    end
    miu(:,:,nt) = miuP;
    P(:,:,:,nt) = PP;
    
    % update continuous state
    meaError = zeros(2,3);
    meaVar = zeros(2,2,3);
    for ns = 1:3
        d = sqrt((miu(1,ns,nt)-xL1)^2+(miu(2,ns,nt)-xL2)^2);
        H = [(miu(1,ns,nt)-xL1)/d,(miu(2,ns,nt)-xL2)/d,0;
            -(miu(2,ns,nt)-xL2)/d^2,(miu(1,ns,nt)-xL1)/d^2,0];
        Hx = [d;atan2(miu(2,ns,nt)-xL2,miu(1,ns,nt)-xL1)];
        meaError(:,ns) = xMea(:,nt)-Hx;
        meaError(2,ns) = wrapToPi(meaError(2,ns));
        meaVar(:,:,ns) = H*P(:,:,ns,nt)*H'+R;
        
        K = P(:,:,ns,nt)*H'*meaVar(:,:,ns)^-1;
        miu(:,ns,nt) = miu(:,ns,nt)+K*meaError(:,ns);
        miu(3,ns,nt) = wrapToPi(miu(3,ns,nt));
        P(:,:,ns,nt) = (eye(3)-K*H)*P(:,:,ns,nt);
    end
    
    % update discrete state
    for ns = 1:3
        ps(ns,nt) = ps(ns,nt)*det(meaVar(:,:,ns))^(-1/2)*...
            exp(-1/2*meaError(:,ns)'*meaVar(:,:,ns)^(-1)*meaError(:,ns));
    end
    ps(:,nt) = ps(:,nt)/sum(ps(:,nt));
    
    % estimate
    xEst(1:2,nt) = miu(1:2,:,nt)*ps(:,nt);
    xEst(3,nt) = atan2(sin(miu(3,:,nt))*ps(:,nt),cos(miu(3,:,nt))*ps(:,nt));
    [~,xEst(4,nt)] = max(ps(:,nt));
    
    tIte(nt-1) = toc(timerIte);
end

tTot = toc(timerTot);

xEst = xEst';

rmpath('..','..\..\lib');

end


function [ R ] = getR( k )

Ns = 1000000;

theta = vmrnd(0,k,Ns);
R = mean(theta.*theta);

end


function [ kernel ] = getKernel( x1, x2, x3, xo1, xo2 )

kernel{1}(:,:,:,1) = wrapToPi(atan2(xo2(1)-x2,xo1(1)-x1')-reshape(x3,1,1,[]))<0 &...
    wrapToPi(atan2(xo2(1)-x2,xo1(1)-x1')-reshape(x3,1,1,[]))>=-pi;
kernel{1}(:,:,:,2) = wrapToPi(atan2(xo2(2)-x2,xo1(2)-x1')-reshape(x3,1,1,[]))<0 &...
    wrapToPi(atan2(xo2(1)-x2,xo1(1)-x1')-reshape(x3,1,1,[]))>=-pi;
kernel{1}(:,:,:,3) = wrapToPi(atan2(xo2(3)-x2,xo1(3)-x1')-reshape(x3,1,1,[]))<0 &...
    wrapToPi(atan2(xo2(1)-x2,xo1(1)-x1')-reshape(x3,1,1,[]))>=-pi;
kernel{2}(:,:,:,1) = wrapToPi(atan2(xo2(1)-x2,xo1(1)-x1')-reshape(x3,1,1,[]))>=0 &...
    wrapToPi(atan2(xo2(1)-x2,xo1(1)-x1')-reshape(x3,1,1,[]))<pi;
kernel{2}(:,:,:,2) = wrapToPi(atan2(xo2(2)-x2,xo1(2)-x1')-reshape(x3,1,1,[]))>=0 &...
    wrapToPi(atan2(xo2(1)-x2,xo1(1)-x1')-reshape(x3,1,1,[]))<pi;
kernel{2}(:,:,:,3) = wrapToPi(atan2(xo2(3)-x2,xo1(3)-x1')-reshape(x3,1,1,[]))>=0 &...
    wrapToPi(atan2(xo2(1)-x2,xo1(1)-x1')-reshape(x3,1,1,[]))<pi;

kernel{1} = double(kernel{1});
kernel{2} = double(kernel{2});

end


function [ PI ] = getPI( miu, P, lamda, kernel, dt, x1, x2, x3 )

fx = zeros(length(x1),length(x2),length(x3),3);
for n1 = 1:length(x1)
    for n2 = 1:length(x2)
        for n3 = 1:length(x3)
            x = [x1(n1);x2(n2);x3(n3)]-miu(:,1);
            fx(n1,n2,n3,1) = 1/sqrt((2*pi)^3*det(P(:,:,1)))*exp(-1/2*x'*P(:,:,1)^-1*x);
            x = [x1(n1);x2(n2);x3(n3)]-miu(:,2);
            fx(n1,n2,n3,2) = 1/sqrt((2*pi)^3*det(P(:,:,2)))*exp(-1/2*x'*P(:,:,2)^-1*x);
            x = [x1(n1);x2(n2);x3(n3)]-miu(:,3);
            fx(n1,n2,n3,3) = 1/sqrt((2*pi)^3*det(P(:,:,3)))*exp(-1/2*x'*P(:,:,3)^-1*x);
        end
    end
end

dx1 = x1(2)-x1(1);
dx2 = x2(2)-x2(1);
dx3 = x3(2)-x3(1);
fx = fx./(sum(sum(sum(fx)))*dx1*dx2*dx3);

PI(1,2) = sum(sum(sum(sum((1-exp(-lamda(:,:,:,1:3)*dt)).*kernel{1}.*fx(:,:,:,1)))))*dx1*dx2*dx3;
PI(1,3) = sum(sum(sum(sum((1-exp(-lamda(:,:,:,1:3)*dt)).*kernel{2}.*fx(:,:,:,1)))))*dx1*dx2*dx3;
PI(2,1) = sum(sum(sum((1-exp(-lamda(:,:,:,4)*dt)).*fx(:,:,:,2))))*dx1*dx2*dx3;
PI(3,1) = sum(sum(sum((1-exp(-lamda(:,:,:,4)*dt)).*fx(:,:,:,3))))*dx1*dx2*dx3;

PI(1,1) = 1-PI(1,2)-PI(1,3);
PI(2,2) = 1-PI(2,1);
PI(2,3) = 0;
PI(3,3) = 1-PI(3,1);
PI(3,2) = 0;

end

