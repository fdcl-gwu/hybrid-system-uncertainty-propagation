clear; close all;
addpath('tests');

% parameters
g = 9.8;                                    % Gravity constant
niu = 0.05;                                 % Air drag coefficient
sigmaNiu = 0.01;                            % standard deviation of air drag coefficient
c = 0.9;                                    % coefficient of restitution
sigmaC = 0.1;                               % standard deviation of coefficient of restitution
epsilonLamda = 0.1;                         % concentration parameter for transition rate
sigmaX1 = 0.2;                              % concentration parameter for position reset
x0 = [1.5;0];                               % initial condition
sigma0 = [0.2^2,0;0,0.5^2];                 % covariance matrix of initial condition

% grid
n1 = 50; n2 = 50;
L1 = 5; L2 = 16;
x1 = linspace(-L1/2,L1/2-L1/n1,n1).'; x1(abs(x1)<1e-10) = 0;
x2 = linspace(-L2/2,L2/2-L2/n2,n2); x2(abs(x2)<1e-10) = 0;
nt = 241;
Lt = 6;
t = linspace(0,Lt,nt);

% initial density function
f0 = @(x)1/(2*pi)/sqrt(det(sigma0))*exp(-1/2*(x-x0)'*sigma0^-1*(x-x0));
fx = zeros(n1,n2,nt);
for m = 1:n1
    for n = 1:n2
        fx(m,n,1) = f0([x1(m);x2(n)]);
    end
end

% initial fft
shift1 = (-1).^((0:n1-1)-floor(n1/2)).';
shift2 = (-1).^((0:n2-1)-floor(n2/2));
y = zeros(n1,n2,nt);
y(:,:,1) = fftshift(fft2(fx(:,:,1)))/n1/n2;
y(:,:,1) = shift1.*shift2.*y(:,:,1);

% fft for f(x)=x and f(x)=x^2
y_x = fftshift(fft(x2))/n2;
y_x = shift2.*y_x;
y_xsqr = fftshift(fft(x2.*abs(x2)))/n2;
y_xsqr = shift2.*y_xsqr;
y_xfour = fftshift(fft(x2.^4))/n2;
y_xfour = shift2.*y_xfour;

% transition kernal and rate
lamda = zeros(n1,n2);
lamda(x1>=0,x2<=0) = repmat(1e3*epsilonLamda*exp(-x1(x1>=0)/epsilonLamda),1,length(find(x2<=0)));
lamda(x1<0,x2<=0) = repmat(lamda(x1==0,x2==0),length(find(x1<0)),length(find(x2<=0)));
kai = zeros(1,n2,n1,n2);
for n_1 = setdiff(1:n2,find(x2==0))
    for m_2 = 1:n1
        for n_2 = 1:n2
            kai(1,n_1,m_2,n_2) = 1/sqrt(2*pi)/sigmaX1*exp(-x1(m_2)^2/2/sigmaX1^2)/(sqrt(2*pi)*sigmaC*abs(x2(n_1)))*...
                exp(-(x2(n_2)+c*x2(n_1))^2/(2*sigmaC^2*x2(n_1)^2));
        end
    end
    kai(1,n_1,:,:) = kai(1,n_1,:,:)/sum(sum(kai(1,n_1,:,:)*L1*L2/n1/n2));
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

% continuous & discrete transition matrix
A = zeros(n1*n2,n1*n2);
for m_1 = 1:n1
    tic;
    for n_1 = 1:n2
        M_1 = m_1-1-floor(n1/2);
        N_1 = n_1-1-floor(n2/2);
        E = L1*L2/(n1*n2)*(exp((-n1/2:n1/2-1).'*M_1/n1).*exp((-n2/2:n2/2-1)*N_1/n2)).^(-2*pi*1i);
        
        for m_2 = 1:n1
            for n_2 = 1:n2
                M_2 = m_2-1-floor(n1/2);
                N_2 = n_2-1-floor(n2/2);
                
                % continuous part
                if m_1 == m_2
                    if m_1 == 1
                        part1 = 0;
                    else
                        part1 = -2*pi*1i*M_1/L1*y_x(wrapDFT(N_1-N_2,n2)+1+floor(n2/2));
                    end
                    
                    if n_1 == 1
                        part2 = 0;
                    else
                        if n_1 == n_2
                            part2 = 2*pi*1i*N_1/L2*g + 2*pi*1i*N_1/L2*niu*y_xsqr(wrapDFT(N_1-N_2,n2)+1+floor(n2/2));
                        else
                            part2 = 2*pi*1i*N_1/L2*niu*y_xsqr(wrapDFT(N_1-N_2,n2)+1+floor(n2/2));
                        end
                    end
                    
                    part3 = -2*sigmaNiu^2*pi^2*N_1^2/L2^2*y_xfour(wrapDFT(N_1-N_2,n2)+1+floor(n2/2));
                else
                    part1 = 0;
                    part2 = 0;
                    part3 = 0;
                end
                
                % discrete part
                part4 = sum(sum(reshape(y_kailamda(wrapDFT(-M_2,n1)+1+floor(n1/2),wrapDFT(-N_2,n2)+1+floor(n2/2),:,:),n1,n2).*E));
                part5 = -y_lamda(wrapDFT(M_1-M_2,n1)+1+floor(n1/2),wrapDFT(N_1-N_2,n2)+1+floor(n2/2));
                
                A((n_1-1)*n1+m_1,(n_2-1)*n1+m_2) = part1+part2+part3+part4+part5;
            end
        end
    end
    toc;
end
expA = expm(A*Lt/(nt-1));

% density reconstruction from Fourier coefficient
fraq1 = ((0:n1-1)'-floor(n1/2))/L1;
fraq2 = ((0:n2-1)-floor(n2/2))/L2;
f_fft = @(x,y)sum(sum(y.*exp(1i*fraq1*2*pi*x(1)).*exp(1i*fraq2*2*pi*x(2))));

% propagation
for i = 2:nt
    y(:,:,i) = reshape(expA*reshape(y(:,:,i-1),[],1),n1,n2);
    for m = 1:n1
        for n = 1:n2
            fx(m,n,i) = real(f_fft([x1(m);x2(n)],y(:,:,i)));
        end
    end
end

% plot
for i = 1:nt
    figure;
    surf(x2,x1,fx(:,:,i));
    view([0,0,1]);
end

% save data
parameter.g = g;
parameter.niu = niu;
parameter.sigmaNiu = sigmaNiu;
parameter.c = c;
parameter.sigmaC = sigmaC;
parameter.sigmaX1 = sigmaX1;
parameter.x0 = x0;
parameter.sigma0 = sigma0;

save(strcat('D:\result-bouncing ball\',sprintf('%i-%i-%i-%i-%i-%i',round(clock)),'.mat'),'parameter','x1','x2','t','y','fx');
