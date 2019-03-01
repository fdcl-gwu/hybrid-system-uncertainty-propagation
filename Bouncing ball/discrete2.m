clear; close all;

% parameters
c = 1;                                      % coefficient of restitution
sigmaC = 0.1;                               % standard deviation of coefficient of restitution
epsilonLamda = 0.1;                         % concentration parameter for transition rate
sigmaX1 = 0.5;                              % concentration parameter for position reset
x0 = [0;-3];                                % initial condition
sigma0 = [0.5^2,0;0,0.5^2];                 % covariance matrix of initial condition

% grid
n1 = 50; n2 = 50;
L1 = 20; L2 = 12;
x1 = linspace(-L1/2,L1/2-L1/n1,n1).'; x1(abs(x1)<1e-10) = 0;
x2 = linspace(-L2/2,L2/2-L2/n2,n2); x2(abs(x2)<1e-10) = 0;
dt = 0.5;

% initial density function
f0 = @(x)1/(2*pi)/sqrt(det(sigma0))*exp(-1/2*(x-x0)'*sigma0^-1*(x-x0));
fx = zeros(n1,n2,2);
for m = 1:n1
    for n = 1:n2
        fx(m,n,1) = f0([x1(m);x2(n)]);
    end
end

% initial fft
shift1 = (-1).^((0:n1-1)-floor(n1/2)).';
shift2 = (-1).^((0:n2-1)-floor(n2/2));
y = zeros(n1,n2,2);
y(:,:,1) = fftshift(fft2(fx(:,:,1)))/n1/n2;
y(:,:,1) = shift1.*shift2.*y(:,:,1);

% transition kernal and rate
lamda = zeros(n1,n2);
lamda(x1>=0,x2<=0) = repmat(1e2*epsilonLamda*exp(-x1(x1>=0)/epsilonLamda),1,length(find(x2<=0)));
lamda(x1<0,x2<=0) = repmat(lamda(x1==0,x2==0),length(find(x1<0)),length(find(x2<=0)));
kai = zeros(1,n2,n1,n2);
for n_1 = setdiff(1:n2,find(x2==0))
    for m_2 = 1:n1
        for n_2 = 1:n2
            kai(1,n_1,m_2,n_2) = 1/sqrt(2*pi)/sigmaX1*exp(-x1(m_2)^2/2/sigmaX1^2)/(sqrt(2*pi)*sigmaC*abs(x2(n_1)))*...
                exp(-(x2(n_2)+c*x2(n_1))^2/(2*sigmaC^2*x2(n_1)^2));
        end
    end
end
kai = repmat(kai,n1,1);

% fft of transition kernal and rate
y_lamda = fftshift(fft2(lamda))/n1/n2;
y_lamda = shift1.*shift2.*y_lamda;
y_kailamda = zeros(n1,n2,n1,n2);
for m_2 = 1:n1
    for n_2 = 1:n2
        y_kailamda(:,:,m_2,n_2) = fftshift(fft2(kai(:,:,m_2,n_2).*lamda))/n1/n2;
        y_kailamda(:,:,m_2,n_2) = shift1.*shift2.*y_kailamda(:,:,m_2,n_2);
    end
end

% propagation
A = zeros(n1*n2,n1*n2);
for m_1 = 1:n1
    for n_1 = 1:n2
        M_1 = m_1-1-floor(n1/2);
        N_1 = n_1-1-floor(n2/2);
        E = L1*L2/(n1*n2)*(exp((-n1/2:n1/2-1).'*M_1/n1).*exp((-n2/2:n2/2-1)*N_1/n2)).^(-2*pi*1i);
        
        for m_2 = 1:n1
            for n_2 = 1:n2
                M_2 = m_2-1-floor(n1/2);
                N_2 = n_2-1-floor(n2/2);
                
                f00 = y_kailamda(wrapDFT(-M_2,n1)+1+floor(n1/2),wrapDFT(-N_2,n2)+1+floor(n2/2),:,:);
                part1 = sum(sum(reshape(f00,n1,n2).*E));
                
                part2 = -y_lamda(wrapDFT(M_1-M_2,n1)+1+floor(n1/2),wrapDFT(N_1-N_2,n2)+1+floor(n2/2));
                
                A((n_1-1)*n1+m_1,(n_2-1)*n1+m_2) = part1+part2;
            end
        end
    end
end

y(:,:,2) = reshape(expm(A*dt)*reshape(y(:,:,1),[],1),n1,n2);

% reconstruct density
fraq1 = ((0:n1-1)'-floor(n1/2))/L1;
fraq2 = ((0:n2-1)-floor(n2/2))/L2;
f_fft = @(x,y)sum(sum(y.*exp(1i*fraq1*2*pi*x(1)).*exp(1i*fraq2*2*pi*x(2))));
for i = 1:n1
    for j = 1:n2
        fx(i,j,2) = real(f_fft([x1(i);x2(j)],y(:,:,2)));
    end
end

% plot
figure;
surf(x2,x1,fx(:,:,1));
view([0,0,1]);

figure;
surf(x2,x1,fx(:,:,2));
view([0,0,1]);

