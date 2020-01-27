function [ x, xMea ] = generateSample( n )

addpath('..','..\..\lib');
close all;

% parameters
Vh = -10;
Vc = 10;
b = 15;
nSample = n;

% grid
Nt = 401;
Lt = 40;
dt = Lt/(Nt-1);

% initial conditions
x1_0 = 220; x2_0 = 0;
s_0 = 1;

% initial measurement
xMea = zeros(2,Nt,nSample);
xMea(1,1,:) = atan2(x2_0,x1_0+9260) + 0.028*sqrt(dt)*randn(1,1,nSample);
xMea(2,1,:) = sqrt((x1_0+9260)^2+x2_0^2) + 40*sqrt(dt)*randn(1,1,nSample);

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
    h2c = poissrnd(lamda(1,:)*dt);
    c2h = poissrnd(lamda(2,:)*dt);
    
    for ns = 1:nSample
        timeh2c = rand(1,h2c(ns));
        timec2h = rand(1,c2h(ns));
        
        if isempty(timeh2c) && isempty(timec2h)
            x(3,nt,ns) = x(3,nt-1,ns);
        elseif isempty(timeh2c)
            x(3,nt,ns) = 1;
        elseif isempty(timec2h)
            x(3,nt,ns) = 2;
        else
            if max(timeh2c)>max(timec2h)
                x(3,nt,ns) = 2;
            else
                x(3,nt,ns) = 1;
            end
        end 
    end
    
    % measurement
    xMea(1,nt,:) = atan2(x(2,nt,:),x(1,nt,:)+9260) + 0.028*sqrt(dt)*randn(1,1,nSample);
    xMea(2,nt,:) = sqrt((x(1,nt,:)+9260).^2+x(2,nt,:).^2) + 40*sqrt(dt)*randn(1,1,nSample);
end

rmpath('..','..\..\lib');

end

