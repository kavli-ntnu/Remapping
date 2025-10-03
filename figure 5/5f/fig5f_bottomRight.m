load simulationVariables_5f_hist

%%
binRange = 0:0.02:1;
counts = histcounts(maxCorrRot,binRange) / sum(~isnan(maxCorrRot));

figure;
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on

counts = histcounts(maxCorrShift,binRange) / sum(~isnan(maxCorrShift));
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on

set(gca,'xtick',0:10:60,'xticklabels',binRange(1:10:end))
