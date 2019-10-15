function [ p ] = getParameter( model )

if model == 1
    % parameters
    p.g = 9.8;
    p.niu = 0.05;
    p.sigma = 0.01;
    p.c = 0.95;
    p.sigmaV = 0.5;
    p.sigmaM = 0.3;
    p.nSample = 1000000;
    
    % initial condition
    p.x0 = [1.5;0];
    p.sigma0 = [0.2^2,0;0,0.5^2];
    
    % grid
    p.N1 = 100;
    p.N2 = 100;
    p.L1 = 5;
    p.L2 = 16;
    p.Nt = 241;
    p.Lt = 6;
end

if model == 2
    % parameters
    p.g = 9.8;
    p.niu = 0.1;
    p.sigma = 0.02;
    p.c = 0.9;
    p.sigmaV = 0.3;
    p.sigmaM = 0.3;
    p.nSample = 1000000;
    
    % initial condition
    p.x0 = [1.5;0];
    p.sigma0 = [0.1^2,0;0,0.2^2];
    
    % grid
    p.N1 = 200;
    p.N2 = 200;
    p.L1 = 5;
    p.L2 = 16;
    p.Nt = 801;
    p.Lt = 20;
end

if model == 3
    % parameters
    p.g = 9.8;
    p.niu = 0.05;
    p.sigma = 0.01;
    p.c = 0.95;
    p.sigmaV = 0.5;
    p.sigmaM = 0.3;
    p.nSample = 10000;
    
    % initial condition
    p.x0 = [1.5;0];
    p.sigma0 = [0.2^2,0;0,0.5^2];
    
    % grid
    p.N1 = 100;
    p.N2 = 100;
    p.L1 = 5;
    p.L2 = 16;
    p.Nt = 241;
    p.Lt = 6;
end

end

