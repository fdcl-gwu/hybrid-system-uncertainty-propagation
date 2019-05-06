function [ x ] = randpdfCar( x1, x2, x3, fx, n )

L1 = x1(end)-x1(1)+x1(2)-x1(1);
L2 = x2(end)-x2(1)+x2(2)-x2(1);
L3 = x3(end)-x3(1)+x3(2)-x3(1);
N1 = length(x1);
N2 = length(x2);
N3 = length(x3);

cdf4 = cumsum(sum(sum(sum(fx,1),2),3),4)*L1/N1*L2/N2*L3/N3;
cdf1 = cumsum(sum(sum(fx,2),3),1)*L1/N1*L2/N2*L3/N3;
cdf2 = cumsum(sum(fx,3),2)*L2/N2*L3/N3;
cdf3 = cumsum(fx,3)*L3/N3;

x = zeros(n,4);
ind = rand(n,4);
for ns = 1:n
    ind4 = ind(ns,4)<cdf4 & ind(ns,4)>=cat(4,0,cdf4(1:end-1));
    x(ns,4) = find(ind4);
    
    ind(ns,1) = ind(ns,1)*cdf1(end,:,:,ind4);
    ind1 = ind(ns,1)<cdf1(:,:,:,ind4) & ind(ns,1)>=cat(1,0,cdf1(1:end-1,:,:,ind4));
    x(ns,1) = x1(ind1);
    
    ind(ns,2) = ind(ns,2)*cdf2(ind1,end,:,ind4);
    ind2 = ind(ns,2)<cdf2(ind1,:,:,ind4) & ind(ns,2)>=cat(2,0,cdf2(ind1,1:end-1,:,ind4));
    x(ns,2) = x2(ind2);
    
    ind(ns,3) = ind(ns,3)*cdf3(ind1,ind2,end,ind4);
    ind3 = ind(ns,3)<cdf3(ind1,ind2,:,ind4) & ind(ns,3)>=cat(3,0,cdf3(ind1,ind2,1:end-1,ind4));
    x(ns,3) = x3(ind3);
end

fxN = zeros(N1,N2,N3,3);
for ns = 1:n
    [~,index1] = min(abs(x(ns,1)-x1));
    [~,index2] = min(abs(x(ns,2)-x2));
    [~,index3] = min(abs(wrapToPi(x(ns,3)-x3)));
    fxN(index1,index2,index3,x(ns,4)) = fxN(index1,index2,index3,x(ns,4))+1;
end
fxN = fxN/n*N1/L1*N2/L2*N3/(2*pi);

end

