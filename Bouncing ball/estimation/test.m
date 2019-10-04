function [] = test()

parpool(12);
for n = 1:12
    rng(n);
end

N = 60;

parfor n = 1:N
    [xTrue,xMea] = generateSample(1);
    xTrue = permute(xTrue,[3,2,1]);
    xMea = xMea';
    
    [xEstS,fxS,tTotS,tIteS,p] = estimateSplitting(xMea,xTrue);
    [xEstMC,fxMC,tTotMC,tIteMC] = estimateMC(xMea,xTrue);
    
    parsave(n,xTrue,xMea,xEstS,fxS,tTotS,tIteS,p,xEstMC,fxMC,tTotMC,tIteMC);
end

end


function [] = parsave(n,xTrue,xMea,xEstS,fxS,tTotS,tIteS,p,xEstMC,fxMC,tTotMC,tIteMC)

save(strcat('D:\result-bouncing ball\10-3-2019\',num2str(n),'.mat'),...
    'xTrue','xMea','fxS','xEstS','tTotS','tIteS','p','fxMC','xEstMC',...
    'tTotMC','tIteMC','-v7.3');

end

