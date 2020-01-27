function [ fx ] = hybrid(  )

close all;
addpath('..','..\..\lib');

% parameters
Vh = -10;
Vc = 10;
b = 15;

% grid
N1 = 120; N2 = 120;
L1 = 600; L2 = 600;
x1 = linspace(-L1/2,L1/2-L1/N1,N1);
x2 = linspace(-L2/2,L2/2-L2/N2,N2);
Nt = 401;
Lt = 40;
t = linspace(0,Lt,Nt);

% initial conditions
x1_0 = 180; x2_0 = 0;
sigma1_0 = 10; sigma2_0 = 10;
s_0 = 1;

% initial distribution
fx = zeros(N1,N2,2,2);
fx(:,:,s_0,1) = 1/(sqrt(2*pi)*sigma1_0)*exp(-0.5*(reshape(x1,[],1)-x1_0).^2/sigma1_0^2) .*...
    (1/sqrt(2*pi)*sigma2_0*exp(-0.5*(reshape(x2,1,[])-x2_0).^2/sigma2_0^2));

% Fourier transform of initial distribution
y = zeros(N1,N2,2);
shift1 = reshape((-1).^((0:N1-1)-N1/2),[],1,1);
shift2 = reshape((-1).^((0:N2-1)-N2/2),1,[],1);
for s = 1:2
    y(:,:,s) = fftshift(fft2(fx(:,:,s,1)))/N1/N2;
    y(:,:,s) = y(:,:,s).*shift1.*shift2;
end

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

% propagation
for nt = 2:Nt
    % continuous
    y = expA.*y;
    
    % reconstruct density
    for s = 1:2
        fx(:,:,s,nt) = ifft2(ifftshift(y(:,:,s)./shift1./shift2*N1*N2),'symmetric');
    end
    
    % renormalize
    fx(:,:,:,nt) = fx(:,:,:,nt)/(sum(sum(sum(fx(:,:,:,nt)))))*L1/N1*L2/N2;
    
    % discrete
    for n1 = 1:N1
        for n2 = 1:N2
            fx(n1,n2,:,nt) = reshape(expB{n1,n2}*reshape(fx(n1,n2,:,nt),2,1,1),1,1,2);
        end
    end
    
    % Fourier transform
    for s = 1:2
        y(:,:,s) = fftshift(fft2(fx(:,:,s,nt)))/N1/N2;
        y(:,:,s) = y(:,:,s).*shift1.*shift2;
    end
end

% plot
for nt = 1:10:Nt
    figure;
    surf(x1,x2,sum(fx(:,:,:,nt),3));
    view([0,0,1]);
end

rmpath('..','..\..\lib');

end

