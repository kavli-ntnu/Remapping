load expList_HP_2d
load placeCorrAll_2d

figure('position',[-743,240,311,420]);
plotSpread(placeCorrAll);
set(findobj(gca,'type','line'),'color','k','marker','o','markersize',4,'linewidth',1)
hold on
plot(1:2,nanmedian(placeCorrAll),'r+','markersize',15,'linewidth',2)
% box on
ylabel('Spatial correlation')
xlim([0.25 2.75])

%%
load placeCorrByExp_2d

placeCorrMean = nan(size(placeCorrByExp,1),2);
for iRow = 1:size(placeCorrByExp,1)
    if ~isempty(placeCorrByExp{iRow,1})
        temp = placeCorrByExp{iRow,1};
        tempMeans = nanmean(temp);
        placeCorrMean(iRow,:) = tempMeans;
    end
end

figure('position',[-1285,67,198,420]);
plotSpread(placeCorrMean)
set(findobj(gca,'type','line'),'color','k','marker','o','markersize',6,'linewidth',1.5)
hold on
plot(1:2,nanmean(placeCorrMean),'r+','markersize',15,'linewidth',2)
set(gca,'xtick',1:2)
xlim([0.5 2.5])
ylim([-0.5 1])
ylabel('Spatial correlation')

%%
load rateDiffAll_2d

figure('position',[-743,240,311,420]);
plotSpread(rateDiffAll)
set(findobj(gca,'type','line'),'color','k','marker','o','markersize',4,'linewidth',1)
hold on
plot(1:2,nanmedian(rateDiffAll),'r+','markersize',15,'linewidth',2)
% box on
ylabel('Rate difference')
xlim([0.25 2.75])
set(gca,'ytick',0:0.25:1)
ylim([-0.02 1])

%%
load rateDiffByExp_2d
rateDiffMean = nan(size(rateDiffByExp,1),2);
for iRow = 1:size(rateDiffByExp,1)
    if ~isempty(rateDiffByExp{iRow,1})
        temp = rateDiffByExp{iRow,1};
        tempMeans = nanmean(temp);
        rateDiffMean(iRow,:) = tempMeans;
    end
end

%
figure('position',[-1285,67,198,420]);
plotSpread(rateDiffMean)
set(findobj(gca,'type','line'),'color','k','marker','o','markersize',6,'linewidth',1.5)
hold on
plot(1:2,nanmean(rateDiffMean),'r+','markersize',15,'linewidth',2)
set(gca,'xtick',1:2)
xlim([0.5 2.5])
ylim([0 1])
set(gca,'ytick',0:0.25:1)
ylabel('Rate difference')
