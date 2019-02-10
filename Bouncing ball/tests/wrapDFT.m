function [ x ] = wrapDFT( x, n )

while x > ceil(n/2)-1
    x = x-n;
end

while x < -floor(n/2)
    x = x+n;
end

end

