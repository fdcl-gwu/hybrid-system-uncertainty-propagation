function [ prob, fx ] = disc(  )
close all;
addpath('..');

% grid
N1 = 120; N2 = 120;
L1 = 600; L2 = 600;
x1 = linspace(-L1/2,L1/2-L1/N1,N1);
x2 = linspace(-L2/2,L2/2-L2/N2,N2);
dt = 0.025;

% initial conditions
x1_0 = -100; x2_0 = 0;
sigma1_0 = 10; sigma2_0 = 10;
s_0 = 1;

% initial distribution
fx = zeros(N1,N2,2,2);
fx(:,:,s_0,1) = 1/(sqrt(2*pi)*sigma1_0)*exp(-0.5*(reshape(x1,[],1)-x1_0).^2/sigma1_0^2) .*...
    (1/sqrt(2*pi)*sigma2_0*exp(-0.5*(reshape(x2,1,[])-x2_0).^2/sigma2_0^2));

% rate and kernel functions
lamda = zeros(N1,N2,2);
for n1 = 1:N1
    for n2 = 1:N2
        if x1(n1) <= -80
            lamda(n1,n2,1) = 15;
            lamda(n1,n2,2) = 1;
        elseif x1(n1) <= 80
            lamda(n1,n2,1) = -7/80*(x1(n1)-80)+1;
            lamda(n1,n2,2) = 7/80*(x1(n1)+80)+1;
        else
            lamda(n1,n2,1) = 1;
            lamda(n1,n2,2) = 15;
        end
    end
end

expB = cell(N1,N2);
for n1 = 1:N1
    for n2 = 1:N2
        B = [-lamda(n1,n2,1), lamda(n1,n2,2);
            lamda(n1,n2,1), -lamda(n1,n2,2)];
        expB{n1,n2} = expm(B*dt);
    end
end

% propagation
for n1 = 1:N1
    for n2 = 1:N2
        fx(n1,n2,:,2) = reshape(expB{n1,n2}*reshape(fx(n1,n2,:,1),2,1,1),1,1,2);
    end
end

prob(1) = sum(sum(fx(:,:,1,2),1),2)*L1/N1*L2/N2;
prob(2) = sum(sum(fx(:,:,2,2),1),2)*L1/N1*L2/N2;

rmpath('..');

end

