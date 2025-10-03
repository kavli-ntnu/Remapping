load fig3de_sessionPairs
load expList_HP_3de
load grpPVstore_ed10c
groupIDs = [1 3 2];

load colorMat
colorMat(4,:) = [0.5 0.5 0.5];

figure('position',[-1719,38,1663,297]);
for iGrp = 1:4
    subplot(1,4,iGrp)
    
    currentGrp = grpPVstore{1,iGrp};
    for i = 1:size(currentGrp,1)
        C = cdfplot(currentGrp{i,1});
        hold on
        set(C,'color',colorMat(iGrp,:))      
    end

    grid off
    title ''
    xlabel('PV correlation')
    ylabel('Frequency')
end

%%
figure;
for iGrp = 1:4
    grpVals = cell2mat(grpPVstore{1,iGrp});
    
    [f,x,flo,fup] = ecdf(grpVals);
    shadedplot(x',flo',fup',colorMat(iGrp,:),colorMat(iGrp,:));
    C = cdfplot(grpVals);
    set(C,'color','k','linewidth',0.5)
    hold on
    
    grid off
    title ''
    xlabel('PV correlation')
    ylabel('Frequency')
end
