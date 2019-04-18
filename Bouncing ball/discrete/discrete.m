function [fx] = discrete()

close all;
addpath('..');

p = getParameter(1);
% parameters
c = p.c;                                    % coefficient of restitution
sigmaV = p.sigmaV;                          % standard deviation of coefficient of restitution
x0 = [0;-3];                                % initial condition
sigma0 = [0.3^2,0.04;0.04,0.5^2];           % covariance matrix of initial condition

% grid
N1 = p.N1; N2 = p.N2;
L1 = p.L1; L2 = p.L2;
x1 = linspace(-L1/2,L1/2-L1/N1,N1).'; x1(abs(x1)<1e-10) = 0;
x2 = linspace(-L2/2,L2/2-L2/N2,N2); x2(abs(x2)<1e-10) = 0;
dt = 1/40;

% initial density function
f0 = @(x)1/(2*pi)/sqrt(det(sigma0))*exp(-1/2*(x-x0)'*sigma0^-1*(x-x0));
fx = zeros(N1,N2,2);
for m = 1:N1
    for n = 1:N2
        fx(m,n,1) = f0([x1(m);x2(n)]);
    end
end

% transition kernal and rate
lamda = zeros(N1,N2);
lamda(x1<0,x2<0) = 100;
lamda(x1==0,x2<0) = 30;
kai = zeros(N1,N2,N1,N2);
for n1 = 1:N1
    for n2 = 1:N2
        for m1 = find(abs(x1-abs(x1(n1)))<1e-3 | [false(N1/2,1);n1==1;false(N1/2-1,1)])
            for m2 = 1:N2
                kai(n1,n2,m1,m2) = N1/L1/(sqrt(2*pi)*sigmaV)*exp(-(x2(m2)+c*x2(n2))^2/(2*sigmaV^2));
            end
        end
        kai(n1,n2,:,:) = kai(n1,n2,:,:)/sum(sum(kai(n1,n2,:,:)*L1*L2/N1/N2));
    end
end

% density propagation
A = zeros(N1*N2,N1*N2);
for n1 = 1:N1
    for j1 = 1:N1
        for n2 = 1:N2
            for j2 = 1:N2
                part1 = kai(n1,n2,j1,j2)*lamda(n1,n2)*L1/N1*L2/N2;
                if n1 == j1 && n2 == j2
                    part2 = -lamda(j1,j2);
                else
                    part2 = 0;
                end
                A((j2-1)*N1+j1,(n2-1)*N1+n1) = part1+part2;
            end
        end
    end
end

fx(:,:,2) = reshape(expm(A*dt)*reshape(fx(:,:,1),[],1),N1,N2);

% plot
figure;
surf(x2,x1,fx(:,:,1));
view([0,0,1]);

figure;
surf(x2,x1,fx(:,:,2));
view([0,0,1]);

rmpath('..');

end

