function [ lamda ] = getLamda( x1, x2, xo1, xo2 )

% x1 and x2 must be column vectors
assert(size(x1,2)==1 && size(x2,2)==1, 'x1 and x2 must be column vectors');

% parameters
epsilonIn = 8;
epsilonOut = 8;
cIn = 400;
cOut = 100;
dOut = 0.5;

% lamda
lamda = zeros(size(x1,1),4);
for no = 1:length(xo1)
    distance = sqrt(sum(([x1,x2]-[xo1(no),xo2(no)]).^2,2));
    lamda(:,no) = cIn*exp(-distance*epsilonIn);
    lamda(:,4) = lamda(:,4)+cOut*exp(-(distance-dOut)*epsilonOut);
end
lamda(:,4) = cOut-lamda(:,4);
lamda(lamda(:,4)<0,4) = 0;

end

