function [ p, fx ] = disc(  )
close all;

% parameters
xo1 = 0;
xo2 = -0.5;
d = 0.5;
epsilon = 1e6;

% grid
N1 = 100; N2 = 100;
L1 = 6; L2 = 6;
x1 = linspace(-L1/2,L1/2-L1/N1,N1);
x2 = linspace(-L2/2,L2/2-L2/N2,N2);
N3 = 50;
x3 = linspace(-pi,pi-2*pi/N3,N3);
dt = 0.025;

% initial conditions
x1_0 = 0; x2_0 = -1;
sigma1_0 = 0.2; sigma2_0 = 0.2;
x3_0 = pi/2;
k_0 = 20;
s_0 = 1;

% initial distribution
fx = zeros(N1,N2,N3,3,2);
fx(:,:,:,s_0,1) = 1/(sqrt(2*pi)*sigma1_0)*exp(-0.5*(reshape(x1,[],1,1)-x1_0).^2/sigma1_0^2) .* ...
    (1/(sqrt(2*pi)*sigma2_0)*exp(-0.5*(reshape(x2,1,[],1)-x2_0).^2/sigma2_0^2)) .* ...
    (1/(2*pi*besseli(0,k_0))*exp(k_0*cos(reshape(x3,1,1,[])-x3_0)));

% rate and kernel functions
thetao = atan2(xo2-x2,xo1-x1');
[G1_1,G1_2,G1_3] = ind2sub([N1,N2,N3],find(sqrt((x1'-xo1).^2+(x2-xo2).^2)<d &...
    wrapToPi(thetao-reshape(x3,1,1,[]))>=0 &...
    wrapToPi(thetao-reshape(x3,1,1,[])<pi)));
[G2_1,G2_2,G2_3] = ind2sub([N1,N2,N3],find(sqrt((x1'-xo1).^2+(x2-xo2).^2)<d &...
    wrapToPi(thetao-reshape(x3,1,1,[]))>=-pi &...
    wrapToPi(thetao-reshape(x3,1,1,[])<0)));
[G3_1,G3_2] = find(sqrt((x1'-xo1).^2+(x2-xo2).^2)>=d);

A = cell(3,1);
A{1} = [-epsilon,0,0;0,0,0;epsilon,0,0];
A{2} = [-epsilon,0,0;epsilon,0,0;0,0,0];
A{3} = [0,epsilon,epsilon;0,-epsilon,0;0,0,-epsilon];
expA = cellfun(@(x)expm(x*dt),A,'UniformOutput',false);

% propagation
for i = 1:length(G1_1)
    fx(G1_1(i),G1_2(i),G1_3(i),:,2) = ...
        reshape(expA{1}*reshape(fx(G1_1(i),G1_2(i),G1_3(i),:,1),[],1),1,1,1,[]);
end

for i = 1:length(G2_1)
    fx(G2_1(i),G2_2(i),G2_3(i),:,2) = ...
        reshape(expA{2}*reshape(fx(G2_1(i),G2_2(i),G2_3(i),:,1),[],1),1,1,1,[]);
end

for i = 1:length(G3_1)
    for n3 = 1:N3
        fx(G3_1(i),G3_2(i),n3,:,2) = ...
            reshape(expA{3}*reshape(fx(G3_1(i),G3_2(i),n3,:,1),[],1),1,1,1,[]);
    end
end

p(1) = sum(sum(sum(fx(:,:,:,1,2),1),2),3)*L1/N1*L2/N2*(2*pi)/N3;
p(2) = sum(sum(sum(fx(:,:,:,2,2),1),2),3)*L1/N1*L2/N2*(2*pi)/N3;
p(3) = sum(sum(sum(fx(:,:,:,3,2),1),2),3)*L1/N1*L2/N2*(2*pi)/N3;

end

