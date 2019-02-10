function [ x ] = wrapToN( x, N )

while x > N
    x = x-(2*N+1);
end

while x < -N
    x = x+(2*N+1);
end

end

