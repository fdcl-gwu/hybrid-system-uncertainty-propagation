clear; close all;

fx = cell(3,2);

% adjust these two numbers line according to "cont" and "contMC"
L1 = 6; L2 = 6;

for mode = 1:3
    fxF = cont(mode);
    fxMC = contMC(mode);
    close all;
    
    fx{mode,1} = fxF(:,:,:,end);
    fx{mode,2} = fxMC(:,:,:,end);
end

[N1,N2,N3] = size(fx{1,1});
x1 = linspace(-L1/2,L1/2-L1/N1,N1);
x2 = linspace(-L2/2,L2/2-L2/N2,N2);
x3 = linspace(-pi,pi-2*pi/N3,N3);

for mode = 1:3
    figure; hold on;
    plot(x3,reshape(sum(sum(fx{mode,1}*L1/N1*L2/N2,1),2),[],1,1));
    plot(x3,reshape(sum(sum(fx{mode,2}*L1/N1*L2/N2,1),2),[],1,1));
    
    figure; hold on;
    surf(x1,x2,sum(fx{mode,1}*2*pi/N3,3)','LineStyle','none','FaceColor','b','FaceAlpha',0.5);
    surf(x1,x2,sum(fx{mode,2}*2*pi/N3,3)','LineStyle','none','FaceColor','r','FaceAlpha',0.5);
    view([0,0,1]);
end