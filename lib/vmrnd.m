function [ theta ] = vmrnd( miu, kappa, n )
% Draw random samples from Von Mises distribution. miu, kappa and n must be
% scalers, theta is a n-by-1 vector
% See Mardia and Jupp, Directional Statistics, 2009, p.43
if ~isscalar(miu) || ~isscalar(kappa)
    error('miu and kappa must be scalars');
end
if ~isscalar(n)
    error('n must be a scalar');
end

theta = zeros(n,1);

a = 1+sqrt(1+4*kappa^2);
b = (a-sqrt(2*a))/(2*kappa);
r = (1+b^2)/(2*b);

unfinished = n;
while unfinished > 0
    % step one
    U1 = rand(unfinished,1);
    z = cos(pi*U1);
    f = (1+r*z)./(r+z);
    c = kappa*(r-f);
    
    % step two
    U2 = rand(unfinished,1);
    success = c.*(2-c)-U2 > 0;
    
    % step three
    unsuccess = log(c./U2)+1-c < 0;
    success = success | ~unsuccess;
    
    % step four
    m = nnz(success);
    U3 = rand(m,1);
    theta(n-unfinished+1:n-unfinished+m) = miu+sign(U3-0.5).*acos(f(success));
    unfinished = unfinished-m;
end

end

