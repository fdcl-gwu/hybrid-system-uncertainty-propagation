function [fx] = hybrid()

close all;
addpath('..','..\..\lib');
tic;

p = getParameter(1);
% parameters
g = p.g;                                    % Gravity constant
niu = p.niu;                                % Air drag coefficient
sigmaNiu = p.sigma;                         % standard deviation of air drag coefficient
c = p.c;                                    % coefficient of restitution
sigmaV = p.sigmaV;                          % standard deviation for velocity reset
x0 = p.x0;                                  % initial condition
sigma0 = p.sigma0;                          % covariance matrix of initial condition

% grid
N1 = p.N1; N2 = p.N2;
L1 = p.L1; L2 = p.L2;
x1 = linspace(-L1/2,L1/2-L1/N1,N1).'; x1(abs(x1)<1e-10) = 0;
x2 = linspace(-L2/2,L2/2-L2/N2,N2); x2(abs(x2)<1e-10) = 0;
Nt = p.Nt;
Lt = p.Lt;
t = linspace(0,Lt,Nt);

% initial density function
f0 = @(x)1/(2*pi)/sqrt(det(sigma0))*exp(-1/2*(x-x0)'*sigma0^-1*(x-x0));
fx = zeros(N1,N2,Nt);
for n1 = 1:N1
    for n2 = 1:N2
        fx(n1,n2,1) = f0([x1(n1);x2(n2)]);
    end
end

% initial fft
shift1 = (-1).^((0:N1-1)-floor(N1/2)).';
shift2 = (-1).^((0:N2-1)-floor(N2/2));
y = zeros(N1,N2,Nt);
y(:,:,1) = fftshift(fft2(fx(:,:,1)))/N1/N2;
y(:,:,1) = shift1.*shift2.*y(:,:,1);

% fft for f(x)=x and f(x)=x^2
y_x = fftshift(fft(x2))/N2;
y_x = shift2.*y_x;
y_xsqr = fftshift(fft(x2.*abs(x2)))/N2;
y_xsqr = shift2.*y_xsqr;
y_xfour = fftshift(fft(x2.^4))/N2;
y_xfour = shift2.*y_xfour;

% A for continuous transition
ACont = zeros(N2,N2,N1);
for n1 = 1:N1
    for n2 = 1:N2
        for m2 = 1:N2
            n2Minusm2 = wrapDFT(n2-m2,N2)+1+N2/2;
            
            if n1 == 1
                part1 = 0;
            else
                part1 = -2*pi*1i*(n1-1-N1/2)/L1*y_x(m2);
            end
            
            if n2 == 2
                part2 = 0;
            else
                if m2-1-N2/2 == 0
                    part2 = 2*pi*1i*(n2-1-N2/2)/L2*g + 2*pi*1i*(n2-1-N2/2)/L2*niu*y_xsqr(m2);
                else
                    part2 = 2*pi*1i*(n2-1-N2/2)/L2*niu*y_xsqr(m2);
                end
            end
            
            part3 = -2*sigmaNiu^2*pi^2*(n2-1-N2/2)^2/L2^2*y_xfour(m2);
            
            ACont(n2,n2Minusm2,n1) = part1 + part2 + part3;
        end
    end
end

expACont = zeros(N2,N2,N1);
for n1 = 1:N1
    expACont(:,:,n1) = expm(ACont(:,:,n1)*Lt/(Nt-1));
end
clear ACont;

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

% A for discrete transition
ADist = zeros(N1*N2,N1*N2);
for n1 = 1:N1
    for m1 = 1:N1
        for n2 = 1:N2
            for m2 = 1:N2
                part1 = kai(n1,n2,m1,m2)*lamda(n1,n2)*L1/N1*L2/N2;
                if n1 == m1 && n2 == m2
                    part2 = -lamda(m1,m2);
                else
                    part2 = 0;
                end
                ADist((m2-1)*N1+m1,(n2-1)*N1+n1) = part1+part2;
            end
        end
    end
end
expADist = expm(ADist*Lt/(Nt-1));

% propagation
for nt = 2:Nt
    % continuous part
    % Fourier coefficient propagation
    for n1 = 1:N1
        y(n1,:,nt) = (expACont(:,:,n1)*y(n1,:,nt-1).').';
    end
    
    % reconstruct density
    fx(:,:,nt) = ifft2(ifftshift(y(:,:,nt)./shift1./shift2*N1*n2),'symmetric');
    
    temp = fx(:,:,nt);
    temp(fx(:,:,nt)<3e-3) = 0;
    fx(:,:,nt) = temp;
    
    % discrete part
    % density propagation
    fx(:,:,nt) = reshape(expADist*reshape(fx(:,:,nt),[],1),N1,N2);
    fx(:,:,nt) = fx(:,:,nt)/(sum(sum(fx(:,:,nt)*L1*L2/N1/N2)));
    
    % fft again
    y(:,:,nt) = fftshift(fft2(fx(:,:,nt)))/N1/N2;
    y(:,:,nt) = shift1.*shift2.*y(:,:,nt);
    
    fprintf(strcat(num2str(nt),'th iteration finished\n'));
end

simulT = toc;

% plot
for nt = 1:Nt
    figure;
    surf(x2,x1,fx(:,:,nt));
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

rmpath('..','..\..\lib');

end
