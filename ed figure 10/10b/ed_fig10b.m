load fig3de_sessionPairs
load expList_HP_3de
load withinSessionStability_ed10b
groupIDs = [1 3 2];

load colorMat

%%
figure('position',[-820,234,711,282]);
for iGrp = 1:size(groupIDs,2)
    subplot(1,3,iGrp)
    
    grpIdx = groupIDs(iGrp);
    inds = find(idx == grpIdx);
    
    ccMean = ccMeanStore(inds,:);
    
    for iRow = 1:size(ccMean,1)
        plot(1:2,ccMean(iRow,:),'-','color',colorMat(iGrp,:),'linewidth',1.5)
        hold on
        plot(1:2,ccMean(iRow,:),'o','color',colorMat(iGrp,:),'linewidth',1.5)
    end
    
    hold on
    plot(1:2,nanmean(ccMean),'k+','markersize',15,'linewidth',2)
    ylim([0.4 1])
    set(gca,'xtick',1:2,'xticklabels',{'Session 1','Session 2'})
    xlim([0.75 2.25])
    groupStore{1,iGrp} = ccMean;
    ylabel('Spatial correlation')
end