function [ prob, fx ] = disc(  )
close all;
addpath('..');

p = getParameter(1);
% parameters
xo1 = p.xo1;
xo2 = p.xo2;
No = length(xo1);

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

% initial distribution
fx = zeros(N1,N2,N3,3,2);
fx(:,:,:,s_0,1) = 1/(sqrt(2*pi)*sigma1_0)*exp(-0.5*(reshape(x1,[],1,1)-x1_0).^2/sigma1_0^2) .* ...
    (1/(sqrt(2*pi)*sigma2_0)*exp(-0.5*(reshape(x2,1,[],1)-x2_0).^2/sigma2_0^2)) .* ...
    (1/(2*pi*besseli(0,k_0))*exp(k_0*cos(reshape(x3,1,1,[])-x3_0)));

% theta
theta = zeros(N1,N2,No);
for no = 1:No
    theta(:,:,no) = atan2(xo2(no)-x2,xo1(no)-x1');
end

% rate and kernel functions
lamdaIn = zeros(N1,N2,3);
lamdaOut = zeros(N1,N2);
for n2 = 1:N2
    lamda = getLamda(x1',ones(N2,1)*x2(n2),xo1,xo2);
    lamdaIn(:,n2,:) = reshape(lamda(:,1:3),N1,1,No);
    lamdaOut(:,n2) = lamda(:,4);
end

expA = cell(N1,N2,N3);
parfor n1 = 1:N1
    for n2 = 1:N2
        for n3 = 1:N3
            t1 = wrapToPi(theta(n1,n2,1)-x3(n3))<0 && wrapToPi(theta(n1,n2,1)-x3(n3))>=-pi;
            t2 = wrapToPi(theta(n1,n2,2)-x3(n3))<0 && wrapToPi(theta(n1,n2,2)-x3(n3))>=-pi;
            t3 = wrapToPi(theta(n1,n2,3)-x3(n3))<0 && wrapToPi(theta(n1,n2,3)-x3(n3))>=-pi;
            A = [-sum(lamdaIn(n1,n2,:),3),lamdaOut(n1,n2),lamdaOut(n1,n2)
                 sum(reshape(lamdaIn(n1,n2,:),1,[]).*[t1,t2,t3]),-lamdaOut(n1,n2),0
                 sum(reshape(lamdaIn(n1,n2,:),1,[]).*[~t1,~t2,~t3]),0,-lamdaOut(n1,n2)];
            expA{n1,n2,n3} = expm(A*dt);
        end
    end
end

% propagation
temp = zeros(N1,N2,N3,3);
parfor n1 = 1:N1
    for n2 = 1:N2
        for n3 = 1:N3
            temp(n1,n2,n3,:) = ...
                reshape(expA{n1,n2,n3}*reshape(fx(n1,n2,n3,:,1),[],1),1,1,1,[]);
        end
    end
end
fx(:,:,:,:,2) = temp;

prob(1) = sum(sum(sum(fx(:,:,:,1,2),1),2),3)*L1/N1*L2/N2*(2*pi)/N3;
prob(2) = sum(sum(sum(fx(:,:,:,2,2),1),2),3)*L1/N1*L2/N2*(2*pi)/N3;
prob(3) = sum(sum(sum(fx(:,:,:,3,2),1),2),3)*L1/N1*L2/N2*(2*pi)/N3;

rmpath('..');

end

