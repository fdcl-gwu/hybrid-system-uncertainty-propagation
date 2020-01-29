function [ x, xMea ] = generateSample( n )

addpath('..','..\..\lib');
close all;

% parameters
Vh = -10;
Vc = 10;
b = 15;
nSample = n;

% grid
Nt = 40001;
Lt = 40;
dt = Lt/(Nt-1);

% initial conditions
x1_0 = 220; x2_0 = -100;
s_0 = 1;

% initial measurement
xMea = zeros(2,Nt,nSample);
xMea(1,1,:) = atan2(x2_0,x1_0+9260) + 0.028*sqrt(0.1)*randn(1,1,nSample);
xMea(2,1,:) = sqrt((x1_0+9260)^2+x2_0^2) + 40*sqrt(0.1)*randn(1,1,nSample);

% pre-allocate memory
x = zeros(3,Nt,nSample);
x(:,1,:) = repmat([x1_0;x2_0;s_0],1,1,nSample);

% propagate samples
for nt = 2:Nt
    % continuous
    Bt = randn(2,1,nSample);
    for ns = 1:nSample
        x(1:2,nt,ns) = x(1:2,nt-1,ns) + [Vh;Vc*(x(3,nt-1,ns)==2)]*dt ...
            + Bt(:,1,nSample)*b*sqrt(dt);
    end
    
    % rate function
    lamda = zeros(2,nSample);
    for ns = 1:nSample
        if x(1,nt-1,ns) <= -80
            lamda(1,ns) = 15;
            lamda(2,ns) = 1;
        elseif x(1,nt-1,ns) <= 80
            lamda(1,ns) = -7/80*(x(1,nt-1,ns)-80)+1;
            lamda(2,ns) = 7/80*(x(1,nt-1,ns)+80)+1;
        else
            lamda(1,ns) = 1;
            lamda(2,ns) = 15;
        end
    end
    
    % discrete
    x(3,nt,:) = x(3,nt-1,:);
    for ns = 1:nSample
        if x(3,nt-1,ns)==1
            if rand > exp(-lamda(1,ns)*dt)
                x(3,nt,ns) = 2;
            end
        else
            if rand > exp(-lamda(2,ns)*dt)
                x(3,nt,ns) = 1;
            end
        end
    end
    
    % measurement
    xMea(1,nt,:) = atan2(x(2,nt,:),x(1,nt,:)+9260) + 0.028*sqrt(0.1)*randn(1,1,nSample);
    xMea(2,nt,:) = sqrt((x(1,nt,:)+9260).^2+x(2,nt,:).^2) + 40*sqrt(0.1)*randn(1,1,nSample);
end

% down sampling
x = x(:,1:100:40001);
xMea = xMea(:,1:100:40001);

rmpath('..','..\..\lib');

end

