load spacingDepth_ed2a
numMods = 3;

for iExp = 1:2
    
    for iMod = 1:3
        modDepth{1,iMod} = modDepthStore{iExp,iMod};
        modSpacing{1,iMod} = mainSpacingStore{iExp,iMod};
    end
    
    load colorMap
    figure;
    for iMod = 1:numMods
        scatter(modDepth{1,iMod},modSpacing{1,iMod},35,colorMap(iMod,:))
        hold on
        box on
    end
    
    figure;
    cnt = 1;
    
    bw = 6;
    
    [f1,x1] = ksdensity(modSpacing{1,1},'bandwidth',bw);
    [pks,loc,width,prom] = findpeaks(f1,x1,'minPeakHeight',1)
    findpeaks(f1,x1,'minPeakHeight',1)
    hold on
    
    [f2,x2] = ksdensity(modSpacing{1,2},'bandwidth',bw);
    [pks,loc,width,prom] = findpeaks(f2,x2,'minPeakHeight',1)
    findpeaks(f2,x2,'minPeakHeight',1)
    
    [f3,x3] = ksdensity(modSpacing{1,3},'bandwidth',bw);
    [pks,loc,width,prom] = findpeaks(f3,x3,'minPeakHeight',1)
    findpeaks(f3,x3,'minPeakHeight',1)
    
    if iExp == 1
        xlim([35 200])
    else
        xlim([35 150])
    end
    ylim([-0.001 0.06])
    set(gca,'ytick',0:0.02:1)
    box on
    
    cnt = cnt + 1;
    
end
