function [fx] = discrete()

close all;

% parameters
c = 0.95;                                   % coefficient of restitution
sigmaV = 0.5;                               % standard deviation of coefficient of restitution
x0 = [0;-3];                                % initial condition
sigma0 = [0.3^2,0.04;0.04,0.5^2];           % covariance matrix of initial condition

% grid
n1 = 100; n2 = 50;
L1 = 5; L2 = 16;
x1 = linspace(-L1/2,L1/2-L1/n1,n1).'; x1(abs(x1)<1e-10) = 0;
x2 = linspace(-L2/2,L2/2-L2/n2,n2); x2(abs(x2)<1e-10) = 0;
dt = 1/40;

% initial density function
f0 = @(x)1/(2*pi)/sqrt(det(sigma0))*exp(-1/2*(x-x0)'*sigma0^-1*(x-x0));
fx = zeros(n1,n2,2);
for m = 1:n1
    for n = 1:n2
        fx(m,n,1) = f0([x1(m);x2(n)]);
    end
end

% transition kernal and rate
lamda = zeros(n1,n2);
lamda(x1<0,x2<=0) = 500;
lamda(x1==0,x2<=0) = 30;
kai = zeros(n1,n2,n1,n2);
for m_1 = 1:n1
    for n_1 = 1:n2
        for m_2 = find(abs(x1-abs(x1(m_1)))<1e-3 | [false(n1/2,1);m_1==1;false(n1/2-1,1)])
            for n_2 = 1:n2
                if x2(n_1) < 0
                    kai(m_1,n_1,m_2,n_2) = n1/L1/(sqrt(2*pi)*sigmaV)*exp(-(x2(n_2)+c*x2(n_1))^2/(2*sigmaV^2));
                else
                    if n_1 == n_2
                        kai(m_1,n_1,m_2,n_2) = n1/L1*n2/L2;
                    end
                end
            end
        end
        kai(m_1,n_1,:,:) = kai(m_1,n_1,:,:)/sum(sum(kai(m_1,n_1,:,:)*L1*L2/n1/n2));
    end
end

% density propagation
A = zeros(n1*n2,n1*n2);
for m_1 = 1:n1
    for m_2 = 1:n1
        for n_1 = 1:n2
            for n_2 = 1:n2
                part1 = kai(m_1,n_1,m_2,n_2)*lamda(m_1,n_1)*L1/n1*L2/n2;
                if m_1 == m_2 && n_1 == n_2
                    part2 = -lamda(m_2,n_2);
                else
                    part2 = 0;
                end
                A((n_2-1)*n1+m_2,(n_1-1)*n1+m_1) = part1+part2;
            end
        end
    end
end

fx(:,:,2) = reshape(expm(A*dt)*reshape(fx(:,:,1),[],1),n1,n2);

% plot
figure;
surf(x2,x1,fx(:,:,1));
view([0,0,1]);

figure;
surf(x2,x1,fx(:,:,2));
view([0,0,1]);

end

