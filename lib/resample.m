function [ xnew, wnew ] = resample( x, w, Ns )

cdf = cumsum(w);

u0 = rand/Ns;
u = u0+(0:(Ns-1))/Ns;

xnew = zeros(Ns,size(x,2));
ind = 1;
for ns = 1:Ns
    while (u(ns)>cdf(ind))
        ind = ind+1;
    end
    xnew(ns,:) = x(ind,:);
end
wnew = ones(Ns,1)/Ns;

end

