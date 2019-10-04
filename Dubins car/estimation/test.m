function [] = test()

parpool(16);
parfor n = 1:16
    rng(n);
end

N = 60;

parfor n = 1:N
    [xTrue,xMea] = generateSample(1);
    xTrue = permute(xTrue,[3,2,1]);
    xMea = permute(xMea,[3,2,1]);
    
    [xEstS,fxS,tTotS,tIteS,p] = estimateSplitting(xMea);
    [xEstMC,fxMC,tTotMC,tIteMC] = estimateMC(xMea);
    [xEstIMM,miuIMM,PIMM,tTotIMM,tIteIMM] = estimateIMM(xMea');
    
    parsave(n,xTrue,xMea,xEstS,fxS,tTotS,tIteS,p,xEstMC,fxMC,...
        tTotMC,tIteMC,xEstIMM,miuIMM,PIMM,tTotIMM,tIteIMM);
end

end


function [] = parsave(n,xTrue,xMea,xEstS,fxS,tTotS,tIteS,p,xEstMC,fxMC,...
    tTotMC,tIteMC,xEstIMM,miuIMM,PIMM,tTotIMM,tIteIMM)

save(strcat('D:\result-dubins car\10-2-2019\',num2str(n),'.mat'),...
    'xTrue','xMea','fxS','xEstS','tTotS','tIteS','p','fxMC','xEstMC',...
    'tTotMC','tIteMC','xEstIMM','miuIMM','PIMM','tTotIMM','tIteIMM','-v7.3');

end