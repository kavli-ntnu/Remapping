load expList_26718
load rotCCstore_26718
numMods = size(rotCCstore,4)
        
%%
for iMod = 1:3
    figure;
    cnt = 1;

    szmap = 40;
    nDeg = 0:3:359;
    selectedDeg = [0:9:99 102:3:135 144:9:351];
    indices = (selectedDeg/3) + 1;
    
    maxCorr = mxStore(iMod,mindStore(iMod));
    
    for ind = 1:size(indices,2)
        subplot(4,size(indices,2)/4,cnt)
        colorMapBRK(rotCCstore(:,:,indices(ind),iMod),'cutoffs',[-inf maxCorr]);
        hold on
        plot([szmap szmap],[0 szmap*2],'k-','linew',1)
        plot([0 szmap*2],[szmap szmap],'k-','linew',1)
    
        title(sprintf('%.2f',mxStore(iMod,indices(ind))))
        cnt = cnt + 1;
    end
end
