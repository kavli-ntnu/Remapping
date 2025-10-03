load cellPairRotationAA_ed11b

binRange = 0:3:60;
[counts,edges] = histcounts(finalWithinCol,binRange);
normCountsW = counts / sum(~isnan(finalWithinCol));
binCenters = edges(1:end-1) + diff(edges)/2;

figure('position',[346,123,793,255]);
subplot(121)
plot(binCenters,normCountsW,'-o')
xlabel('Rotation difference (°)')
ylabel('Frequency')
title(sprintf('Median = %.1f deg',nanmedian(finalWithinCol)))
box on
xlim([-1 60])
set(gca,'ytick',0:0.1:0.5)

[counts,edges] = histcounts(finalBetweenCol,binRange);
normCountsB = counts / sum(~isnan(finalBetweenCol));
binCenters = edges(1:end-1) + diff(edges)/2;

subplot(122)
plot(binCenters,normCountsB,'-o')
xlabel('Rotation difference (°)')
ylabel('Frequency')
title(sprintf('Median = %.1f deg',nanmedian(finalBetweenCol)))
box on

set(gca,'ytick',0:0.1:0.5)
% ylim([0 0.3])
xlim([-1 60])

%%
load cellPairRotationAB_ed11b

binRange = 0:3:60;
[counts,edges] = histcounts(finalWithinCol,binRange);
normCountsW = counts / sum(~isnan(finalWithinCol));
binCenters = edges(1:end-1) + diff(edges)/2;

figure('position',[346,123,793,255]);
subplot(121)
plot(binCenters,normCountsW,'-o')
xlabel('Rotation difference (°)')
ylabel('Frequency')
title(sprintf('Median = %.1f deg',nanmedian(finalWithinCol)))
box on
xlim([-1 60])
ylim([0 0.3])
set(gca,'ytick',0:0.1:0.5)

[counts,edges] = histcounts(finalBetweenCol,binRange);
normCountsB = counts / sum(~isnan(finalBetweenCol));
binCenters = edges(1:end-1) + diff(edges)/2;

subplot(122)
plot(binCenters,normCountsB,'-o')
xlabel('Rotation difference (°)')
ylabel('Frequency')
title(sprintf('Median = %.1f deg',nanmedian(finalBetweenCol)))
box on

set(gca,'ytick',0:0.1:0.5)
ylim([0 0.3])
xlim([-1 60])
