function [ p ] = getParameter( model )

if model == 1
    % parameters
    p.v = 1;
    p.u = [0,2,-2];
    p.sigma = 0.2;
    
    p.xo1 = [0,-1.5,1.5];
    p.xo2 = [0,1,1];
    
    p.xL1 = 0;
    p.xL2 = -3;
    p.sigmaL = 0.5;
    p.kL = 30;
    
    p.nSample = 1000000;

    % grid
    p.N1 = 100;
    p.L1 = 6;
    p.N2 = 100;
    p.L2 = 6;
    p.N3 = 50;
    p.Nt = 161;
    p.Lt = 4;

    % initial conditions
    p.x1_0 = 0; p.x2_0 = -2;
    p.sigma1_0 = 0.2; p.sigma2_0 = 0.2;
    p.x3_0 = pi/2;
    p.k_0 = 20;
    p.s_0 = 1;
end

end

