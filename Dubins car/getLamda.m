function [ lamda ] = getLamda( x1, x2, xo1, xo2 )

% x1 and x2 must be column vectors
assert(size(x1,2)==1 && size(x2,2)==1, 'x1 and x2 must be column vectors');

% parameters
dIn1 = 0.1;
dIn2 = 0.5;
dOut1 = 0.5;
dOut2 = 0.9;
peak = 50;

% lamda
lamda = zeros(size(x1,1),4);
for no = 1:length(xo1)
    distance = sqrt(sum(([x1,x2]-[xo1(no),xo2(no)]).^2,2));
    
    % In
    lamda(distance<dIn1,no) = peak;
    lamda(distance>=dIn2,no) = 0;
    lamda(distance>=dIn1 & distance<dIn2,no) = ...
        peak*sin((dIn2-distance(distance>=dIn1 & distance<dIn2))/(dIn2-dIn1)*pi/2);
    
    % Out
    lamda(distance>dOut2,4) = lamda(distance>dOut2,4)+peak;
    lamda(distance<=dOut1,4) = lamda(distance<dOut1,4);
    lamda(distance>dOut1 & distance<=dOut2,4) = lamda(distance>dOut1 & distance<=dOut2,4) + ...
        peak*sin((distance(distance>dOut1 & distance<=dOut2)-dOut1)/(dOut2-dOut1)*pi/2);
end
lamda(:,4) = lamda(:,4)-peak*2;

end

