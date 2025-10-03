
load expList_HP_2f
expCol = [3 3 1 3 3 3 3 3];

pvMat = [];
pvMedian = nan(size(expList_HP,1),1);
for iExp = 1:size(expList_HP,1)
    load(sprintf('%s',expList_HP{iExp,8}))
    col = expCol(1,iExp);
    
    pvMat = [pvMat; simplePVstoreRaw(:,col,2)];
    pvMedian(iExp,1) = nanmedian(simplePVstoreRaw(:,col,2));
end

%%
expCol = [1 1 0 2 1 1 1 1];

pvStab = [];
pvStabMedian = nan(size(expList_HP,1),1);
for iExp = 1:size(expList_HP,1)
    if iExp == 3
        pvStabMedian(iExp,1) = nan;
    else
        load(sprintf('%s',expList_HP{iExp,8}))
        col = expCol(1,iExp);
        
        pvStab = [pvStab; simplePVstoreRaw(:,col,2)];
        pvStabMedian(iExp,1) = nanmedian(simplePVstoreRaw(:,col,2));
    end
end

%%
figure;
cdfplot(pvStab)
hold on
cdfplot(pvMat)
title ''
grid off
xlim([-0.5 1])
set(gca,'ytick',0:0.2:1)
ylabel('PV correlation')

%%
figure('position',[-1285,67,138,420]);
plotSpread({pvStabMedian pvMedian})
set(findobj(gca,'type','line'),'color','k','marker','o','markersize',5,'linewidth',1)
plot(1:2,[nanmean(pvStabMedian) nanmean(pvMedian)],'r+','markersize',15,'linewidth',2)
set(gca,'xtick',1:2)
xlim([0.5 2.5])
ylabel('PV correlation')



    
    
    
    