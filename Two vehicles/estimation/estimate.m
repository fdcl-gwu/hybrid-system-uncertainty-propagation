function [ xEst, fx, tTot, tIte ] = estimate( xMea )

close all;
addpath('..','..\..\lib');

timerTot = tic;

% parameters
Vh = -10;
Vc = 10;
b = 15;
sigmaL = [0.028;40];

% grid
N1 = 120; N2 = 120;
L1 = 600; L2 = 600;
x1 = linspace(-L1/2,L1/2-L1/N1,N1);
x2 = linspace(-L2/2,L2/2-L2/N2,N2);
Nt = 401;
Lt = 40;
t = linspace(0,Lt,Nt);

% initial conditions
x1_0 = 220; x2_0 = -100;
sigma1_0 = 10; sigma2_0 = 10;
s_0 = 1;

% initial knowledge
fx = zeros(N1,N2,2,Nt);
fx(:,:,s_0,1) = 1/(sqrt(2*pi)*sigma1_0)*exp(-0.5*(reshape(x1,[],1)-x1_0).^2/sigma1_0^2) .*...
    (1/(sqrt(2*pi)*sigma2_0)*exp(-0.5*(reshape(x2,1,[])-x2_0).^2/sigma2_0^2));

% Fourier transform of initial distribution
y = zeros(N1,N2,2);
shift1 = reshape((-1).^((0:N1-1)-N1/2),[],1,1);
shift2 = reshape((-1).^((0:N2-1)-N2/2),1,[],1);
for s = 1:2
    y(:,:,s) = fftshift(fft2(fx(:,:,s,1)))/N1/N2;
    y(:,:,s) = y(:,:,s).*shift1.*shift2;
end

% likelihood function
sigmaL = sigmaL*sqrt(Lt/(Nt-1));
l = @(alpha,d,x1,x2) (1/2*pi/sigmaL(1)/sigmaL(2))*exp(-(atan2(x2,x1)-alpha)^2/2/sigmaL(1)^2)...
    *exp(-(sqrt(x1^2+x2^2)-d)^2/2/sigmaL(2)^2);

% coefficients continuous propagation
expA = zeros(N1,N2,2);
c1 = [0,ones(1,N1-1)];
c2 = [0,ones(1,N2-1)];
for s = 1:2
    for n1 = 1:N1
        for n2 = 1:N2
            part1 = -2*pi*1i*(n1-1-N1/2)*c1(n1)*Vh/L1;
            part2 = -2*pi*1i*(n2-1-N2/2)*c2(n2)*(s==2)*Vc/L2;
            part3 = -0.5*b^2*(4*pi^2*(n1-1-N1/2)^2/L1^2 + 4*pi^2*(n2-1-N2/2)^2/L2^2);
            A = part1+part2+part3;
            expA(n1,n2,s) = exp(A*Lt/(Nt-1));
        end
    end
end

% coefficients discrete propagation
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
        expB{n1,n2} = expm(B*Lt/(Nt-1));
    end
end

% pre-allocate memory
xEst = zeros(Nt,3);
xEst(1,:) = [x1_0,x2_0,s_0];
tIte = zeros(Nt-1,1);

% estimation
for nt = 2:Nt
    timerIte = tic;
    
    % propagation
    y = expA.*y;
    
    for s = 1:2
        fx(:,:,s,nt) = ifft2(ifftshift(y(:,:,s)./shift1./shift2*N1*N2),'symmetric');
    end
    
    for n1 = 1:N1
        for n2 = 1:N2
            fx(n1,n2,:,nt) = reshape(expB{n1,n2}*reshape(fx(n1,n2,:,nt),2,1,1),1,1,2);
        end
    end

    fx(fx<0) = 0;
    
    % measurement update
    lx = zeros(N1,N2,2);
    for n1 = 1:N1
        for n2 = 1:N2
            lx(n1,n2,:,:) = l(xMea(nt,1),xMea(nt,2),x1(n1)+9260,x2(n2));
        end
    end
    
    fx(:,:,:,nt) = fx(:,:,:,nt).*lx;
    fx(:,:,:,nt) = fx(:,:,:,nt)/(sum(sum(sum(fx(:,:,:,nt))))*L1/N1*L2/N2);
    
    for s = 1:2
        y(:,:,s) = fftshift(fft2(fx(:,:,s,nt)))/N1/N2;
        y(:,:,s) = y(:,:,s).*shift1.*shift2;
    end
    
    % estimation
    xEst(nt,1) = sum(x1'.*sum(sum(fx(:,:,:,nt),2),3)*L2/N2)*L1/N1;
    xEst(nt,2) = sum(x2.*sum(sum(fx(:,:,:,nt),1),3)*L1/N1)*L2/N2;
    [~,xEst(nt,3)] = max(sum(sum(fx(:,:,:,nt),1),2));
    
    tIte(nt-1) = toc(timerIte);
end

tTot = toc(timerTot);

% % plot
% for nt = 1:10:Nt
%     figure;
%     surf(x1,x2,sum(fx(:,:,:,nt),3));
%     view([0,0,1]);
% end

rmpath('..','..\..\lib');

end

