function [fx] = hybrid()

close all;
addpath('tests');
tic;

% parameters
g = 9.8;                                    % Gravity constant
niu = 0.05;                                 % Air drag coefficient
sigmaNiu = 0.01;                            % standard deviation of air drag coefficient
c = 0.95;                                   % coefficient of restitution
sigmaV = 0.5;                               % standard deviation for velocity reset
x0 = [1.5;0];                               % initial condition
sigma0 = [0.2^2,0;0,0.5^2];                 % covariance matrix of initial condition

% grid
n1 = 100; n2 = 50;
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

% A for continuous transition
expACont = zeros(n2,n2,n1);
for m = 1:n1
    ACont = zeros(n2,n2);
    M = m-1-floor(n1/2);
    for n = 1:n2
        N = n-1-floor(n2/2);
        for k = 1:n2
            K = k-1-floor(n2/2);
            NMinusK = wrapDFT(N-K,n2);
            nMinusk = NMinusK+1+floor(n2/2);

            if m == 1 && mod(n1,2) == 0
                part1 = 0;
            else
                part1 = -2*pi*1i*M/L1*y_x(k);
            end

            if n == 1 && mod(n2,2) == 0
                part2 = 0;
            else
                if K == 0
                    part2 = 2*pi*1i*N/L2*g + 2*pi*1i*N/L2*niu*y_xsqr(k);
                else
                    part2 = 2*pi*1i*N/L2*niu*y_xsqr(k);
                end
            end

            part3 = -2*sigmaNiu^2*pi^2*N^2/L2^2*y_xfour(k);

            ACont(n,nMinusk) = part1 + part2 + part3;
        end
    end
    
    expACont(:,:,m) = expm(ACont*Lt/(nt-1));
end

% density reconstruction from Fourier coefficient
fraq1 = ((0:n1-1)'-floor(n1/2))/L1;
fraq2 = ((0:n2-1)-floor(n2/2))/L2;
f_fft = @(x,y)sum(sum(y.*exp(1i*fraq1*2*pi*x(1)).*exp(1i*fraq2*2*pi*x(2))));

% transition kernal and rate
lamda = zeros(n1,n2);
lamda(x1<0,x2<0) = 500;
lamda(x1==0,x2<0) = 30;
kai = zeros(n1,n2,n1,n2);
for m_1 = 1:n1
    for n_1 = 1:n2
        for m_2 = find(abs(x1-abs(x1(m_1)))<1e-3 | [false(n1/2,1);m_1==1;false(n1/2-1,1)])
            for n_2 = 1:n2
                kai(m_1,n_1,m_2,n_2) = n1/L1/(sqrt(2*pi)*sigmaV)*exp(-(x2(n_2)+c*x2(n_1))^2/(2*sigmaV^2));
            end
        end
        kai(m_1,n_1,:,:) = kai(m_1,n_1,:,:)/sum(sum(kai(m_1,n_1,:,:)*L1*L2/n1/n2));
    end
end

% A for discrete transition
ADist = zeros(n1*n2,n1*n2);
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
                ADist((n_2-1)*n1+m_2,(n_1-1)*n1+m_1) = part1+part2;
            end
        end
    end
end
expADist = expm(ADist*Lt/(nt-1));

% propagation
for i = 2:nt
    % continuous part
    % Fourier coefficient propagation
    for m = 1:n1
        y(m,:,i) = (expACont(:,:,m)*y(m,:,i-1).').';
    end
    
    % reconstruct density
    for m = 1:n1
        for n = 1:n2
            fx(m,n,i) = real(f_fft([x1(m);x2(n)],y(:,:,i)));
        end
    end
    
    temp = fx(:,:,i);
    temp(fx(:,:,i)<3e-3) = 0;
    fx(:,:,i) = temp;
    
    % discrete part
    % density propagation
    fx(:,:,i) = reshape(expADist*reshape(fx(:,:,i),[],1),n1,n2);
    fx(:,:,i) = fx(:,:,i)/(sum(sum(fx(:,:,i)*L1*L2/n1/n2)));
    
    % fft again
    y(:,:,i) = fftshift(fft2(fx(:,:,i)))/n1/n2;
    y(:,:,i) = shift1.*shift2.*y(:,:,i);
    
    fprintf(strcat(num2str(i),'th iteration finished\n'));
end

simulT = toc;

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
parameter.sigmaV = sigmaV;
parameter.x0 = x0;
parameter.sigma0 = sigma0;

save(strcat('D:\result-bouncing ball\',sprintf('%i-%i-%i-%i-%i-%i',round(clock)),'.mat'),'parameter','x1','x2','t','y','fx','simulT');

end
