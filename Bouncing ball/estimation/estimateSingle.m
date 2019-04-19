function [ fx, xTrue, xEst ] = estimateSingle(  )
close all;
rng('shuffle');
addpath('..','..\..\lib');

p = getParameter(1);
% parameters
g = p.g;                                    % Gravity constant
niu = p.niu;                                % Air drag coefficient
sigmaNiu = p.sigma;                         % standard deviation of air drag coefficient
c = p.c;                                    % coefficient of restitution
sigmaV = p.sigmaV;                          % standard deviation for velocity reset
sigmaM = p.sigmaM;                          % standard deviation of measurement

% grid
N1 = p.N1; N2 = 100;
L1 = p.L1; L2 = p.L2;
x1 = linspace(-L1/2,L1/2-L1/N1,N1).'; x1(abs(x1)<1e-10) = 0;
x2 = linspace(-L2/2,L2/2-L2/N2,N2); x2(abs(x2)<1e-10) = 0;
Nt = p.Nt;
Lt = p.Lt;
t = linspace(0,Lt,Nt);

% initial knowledge
fx = zeros(N1,N2,Nt);
fx(:,:,1) = 1/L1/L2;

% initial fft
shift1 = (-1).^((0:N1-1)-floor(N1/2)).';
shift2 = (-1).^((0:N2-1)-floor(N2/2));
y = zeros(N1,N2,Nt);
y(:,:,1) = fftshift(fft2(fx(:,:,1)))/N1/N2;
y(:,:,1) = shift1.*shift2.*y(:,:,1);

% true state
xTrue = generateSample(1,p);
xTrue = reshape(xTrue,2,Nt)';

% likelihood function
l = @(h,x1) 1/sqrt(2*pi)/sigmaM*exp(-(x1-h).^2/2/sigmaM^2);

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

% estimation
xEst = zeros(Nt,2);
for nt = 2:Nt
    % propagation
    for n1 = 1:N1
        y(n1,:,nt) = (expACont(:,:,n1)*y(n1,:,nt-1).').';
    end

    fx(:,:,nt) = ifftn(ifftshift(y(:,:,nt)./shift1./shift2*N1*n2),'symmetric');
    
    temp = fx(:,:,nt);
    temp(fx(:,:,nt)<3e-3) = 0;
    fx(:,:,nt) = temp;
    
    fx(:,:,nt) = reshape(expADist*reshape(fx(:,:,nt),[],1),N1,N2);
    
    % measurement update
    lx = repmat(l(xTrue(nt,1),x1),1,N2);
    fx(:,:,nt) = fx(:,:,nt).*lx;
    fx(:,:,nt) = fx(:,:,nt)/(sum(sum(fx(:,:,nt)*L1*L2/N1/N2)));
    
    y(:,:,nt) = fftshift(fft2(fx(:,:,nt)))/N1/N2;
    y(:,:,nt) = shift1.*shift2.*y(:,:,nt);
    
    % estimation
    [~,index1] = max(max(fx(:,:,nt),[],2),[],1);
    [~,index2] = max(max(fx(:,:,nt),[],1),[],2);
    xEst(nt,:) = [x1(index1),x2(index2)];
end

% plot
for nt = 1:4:Nt
    figure; hold on;
    surf(x2,x1,fx(:,:,nt));
    scatter3(xTrue(nt,2),xTrue(nt,1),max(max(fx(:,:,nt))),'MarkerFaceColor','r','MarkerEdgeColor','r');
    scatter3(xEst(nt,2),xEst(nt,1),max(max(fx(:,:,nt))),'MarkerFaceColor','b','MarkerEdgeColor','b');
    view([0,0,1]);
end

save(strcat('D:\result-bouncing ball\',sprintf('%i-%i-%i-%i-%i-%i',round(clock)),'-estimate.mat'),'fx','xTrue','xEst');

rmpath('..','..\..\lib');

end

