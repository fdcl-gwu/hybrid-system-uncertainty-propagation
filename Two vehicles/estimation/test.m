function [  ] = test(  )

parpool(12);
parfor n = 1:12
    rng(n);
end

N = 60;

parfor n = 1:N
    [xTrue,xMea] = generateSample(1);
    xTrue = xTrue';
    xMea = xMea';
    
    [xEst,fx,tTot,tIte] = estimate(xMea);
    
    parsave(n,xTrue,xMea,xEst,fx,tTot,tIte);
end

end


function [] = parsave(n,xTrue,xMea,xEst,fx,tTot,tIte)

save(strcat('D:\result-hybrid\two vehicle\1-26-2020-2\',num2str(n),'.mat'),...
    'xTrue','xMea','xEst','fx','tTot','tIte');

end

