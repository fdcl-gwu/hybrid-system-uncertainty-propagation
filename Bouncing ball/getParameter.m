function [ p ] = getParameter( model )

if model == 1
    % parameters
    p.g = 9.8;
    p.niu = 0.05;
    p.sigma = 0.01;
    p.c = 0.95;
    p.sigmaV = 0.5;
    p.nSample = 1000000;
    
    % initial condition
    p.x0 = [1.5;0];
    p.sigma0 = [0.2^2,0;0,0.5^2];
    
    % grid
    p.N1 = 100;
    p.N2 = 50;
    p.L1 = 5;
    p.L2 = 16;
    p.Nt = 241;
    p.Lt = 6;
end

end

