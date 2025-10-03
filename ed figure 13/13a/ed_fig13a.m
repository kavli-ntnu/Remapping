
load remappingMetrics_translationSims

figure('position',[914,161,1855,312]);
subplot(141)
binRange = 0:0.04:1;
counts = histcounts(maxCorrAll,binRange) / sum(~isnan(maxCorrAll));
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on
% xline(find(binRange == 0.20))
set(gca,'xtick',0:10:60,'xticklabels',binRange(1:10:end))
xlabel('PV correlation')
ylabel('Proportion')
title('Translation')

subplot(142)
binRange = 0:0.04:1.2;
counts = histcounts(rValAll,binRange) / sum(~isnan(rValAll));
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on
% xline(find(binRange == 0.20))
set(gca,'xtick',0:10:60,'xticklabels',binRange(1:10:end))
xlabel('Rearrangement score')

subplot(143)
binRange = 0:0.03:0.74;
counts = histcounts(validStoreAll,binRange) / sum(~isnan(validStoreAll));
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on
% xline(find(binRange == 0.20))
set(gca,'xtick',0:10:60,'xticklabels',binRange(1:10:end))
xlabel('Percent turnover')

subplot(144)
binRange = 0:0.04:1;
counts = histcounts(minDist,binRange) / sum(~isnan(minDist));
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on
% xline(find(binRange == 0.20))
set(gca,'xtick',0:10:60,'xticklabels',binRange(1:10:end))
xlabel('Minimum distance')


%% 
load remappingMetrics_rotationSims

figure('position',[914,161,1855,312]);
subplot(141)
binRange = 0:0.04:1;
counts = histcounts(maxCorrAll,binRange) / sum(~isnan(maxCorrAll));
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on
% xline(find(binRange == 0.20))
set(gca,'xtick',0:10:60,'xticklabels',binRange(1:10:end))
xlabel('PV correlation')
ylabel('Proportion')
title('Rotation')

subplot(142)
binRange = 0:0.04:1.2;
counts = histcounts(remappingVal,binRange) / sum(~isnan(remappingVal));
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on
% xline(find(binRange == 0.20))
set(gca,'xtick',0:10:60,'xticklabels',binRange(1:10:end))
xlabel('Rearrangement score')

subplot(143)
binRange = 0:0.03:0.74;
counts = histcounts(validStoreAll,binRange) / sum(~isnan(validStoreAll));
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on
% xline(find(binRange == 0.20))
set(gca,'xtick',0:10:60,'xticklabels',binRange(1:10:end))
xlabel('Percent turnover')

subplot(144)
binRange = 0:3:30;
counts = histcounts(minDiff,binRange) / sum(~isnan(minDiff));
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on
% xline(find(binRange == 0.20))
set(gca,'xtick',0.5:10.5,'xticklabels',(0:10)*3)
xlabel('Minimum angular difference')

%%
load remappingMetrics_shuffleSims

figure('position',[914,161,1855,312]);
subplot(141)
binRange = 0:0.04:1;
counts = histcounts(maxCorr,binRange) / sum(~isnan(maxCorr));
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on
set(gca,'xtick',0:10:60,'xticklabels',binRange(1:10:end))
xlabel('PV correlation')
ylabel('Proportion')
title('Shuffle')

prc = prctile(maxCorr,95);
[idx1,idx2] = min(abs(binRange - prc));
xline(idx2)

subplot(142)
binRange = 0:0.04:1.2;
counts = histcounts(rVal,binRange) / sum(~isnan(rVal));
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on
set(gca,'xtick',0:10:60,'xticklabels',binRange(1:10:end))
xlabel('Rearrangement score')

prc = prctile(rVal,5);
[idx1,idx2] = min(abs(binRange - prc));
xline(idx2)

subplot(143)
binRange = 0:0.03:0.74;
counts = histcounts(validStoreShuf,binRange) / sum(~isnan(validStoreShuf));
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on
set(gca,'xtick',0:10:60,'xticklabels',binRange(1:10:end))
xlabel('Percent turnover')

prc = prctile(validStoreShuf,5);
[idx1,idx2] = min(abs(binRange - prc));
xline(idx2)

load minDistShuf

subplot(144)
binRange = 0:0.04:0.92;
counts = histcounts(minMat,binRange) / sum(~isnan(minMat));
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on
set(gca,'xtick',0:10:60,'xticklabels',binRange(1:10:end))
xlabel('Minimum distance')

prc = prctile(minMat,5);
[idx1,idx2] = min(abs(binRange - prc));
xline(idx2)

%%
load remappingMetrics_simRotShift

figure('position',[914,161,1855,312]);
subplot(141)
binRange = 0:0.04:1;
counts = histcounts(maxCorr,binRange) / sum(~isnan(maxCorr));
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on
set(gca,'xtick',0:10:60,'xticklabels',binRange(1:10:end))
xlabel('PV correlation')
ylabel('Proportion')
title('Combined')


subplot(142)
binRange = 0:0.04:1.2;
counts = histcounts(rVal,binRange) / sum(~isnan(rVal));
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on
set(gca,'xtick',0:10:60,'xticklabels',binRange(1:10:end))
xlabel('Rearrangement score')

subplot(143)
binRange = 0:0.03:0.74;
counts = histcounts(validStoreAll,binRange) / sum(~isnan(validStoreAll));
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on
set(gca,'xtick',0:10:60,'xticklabels',binRange(1:10:end))
xlabel('Percent turnover')

subplot(144)
binRange = 0:0.04:0.92;
counts = histcounts(minDist,binRange) / sum(~isnan(minDist));
bar(1:size(binRange,2)-1,counts,'barwidth',1,'edgecolor','k','facealpha',0.5)
hold on
set(gca,'xtick',0:10:60,'xticklabels',binRange(1:10:end))
xlabel('Minimum distance')
