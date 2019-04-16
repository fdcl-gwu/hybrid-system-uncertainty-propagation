function [ lamda ] = getLamda( x1, x2, xo1, xo2 )

% x1 and x2 must be column vectors
assert(size(x1,2)==1 && size(x2,2)==1, 'x1 and x2 must be column vectors');

% parameters
epsilonIn = 6;
cIn = 400;
epsilonOut = 0.3;
cOut = 100;

% lamda
lamda = [zeros(size(x1,1),3),cOut*ones(size(x1,1),4)];
for no = 1:length(xo1)
    distance = sqrt(sum(([x1,x2]-[xo1(no),xo2(no)]).^2,2));
    lamda(:,no) = cIn*exp(-distance*epsilonIn);
    lamda(:,4) = lamda(:,4).*exp(-epsilonOut./distance);
end

end

