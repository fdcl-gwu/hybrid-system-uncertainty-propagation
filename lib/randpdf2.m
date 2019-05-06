function [ x ] = randpdf2( x1, x2, fx, n )
% draw random samples from a custom pdf, fx must be normalized

L1 = x1(end)-x1(1)+x1(2)-x1(1);
L2 = x2(end)-x2(1)+x2(2)-x2(1);
N1 = length(x1);
N2 = length(x2);

cdf1 = sum(cumsum(fx,1),2)*L1/N1*L2/N2;
cdf2 = cumsum(fx,2)*L2/N2;

x = zeros(n,2);
ind = rand(n,2);
for ns = 1:n
    ind1 = ind(ns,1)<cdf1 & ind(ns,1)>=[0;cdf1(1:end-1)];
    x(ns,1) = x1(ind1);
    ind(ns,2) = ind(ns,2)*cdf2(ind1,end);
    x(ns,2) = x2(ind(ns,2)<cdf2(ind1,:) & ind(ns,2)>=[0,cdf2(ind1,1:end-1)]);
end

end

