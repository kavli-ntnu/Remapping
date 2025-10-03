%%
load spacing_ed2c

numMods = 3;
m = nan(50,3);

colorblindMap
figure('position',[723,214,274,314]);
hold on
for iExp = 1:size(modMeans,1)
    for iMod = 1:numMods
        currentSpacing = mainSpacingStore{iExp,iMod};
        errorbar(iMod,nanmean(currentSpacing),nanstd(currentSpacing) ./ sqrt(sum(~isnan(currentSpacing))),'color',colorMap(iMod,:),'linewidth',1.5)
    end
    hold on
    plot(1:numMods,modMeans(iExp,:),'-','color',[0.5 0.5 0.5],'linewidth',1.5)  
    m(iExp,:) = modMeans(iExp,:);
end
ylim([40 175])
xlim([0.75 3.25])
ylabel('Grid spacing')
set(gca,'xtick',1:4,'xticklabels',{'M1','M2','M3'})
box on

