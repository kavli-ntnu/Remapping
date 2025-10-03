load fieldDist_ed4c

binRange = 0:3:60;
[counts,edges] = histcounts(fieldDistConcatAA,binRange);
normCountsAA = counts / sum(~isnan(fieldDistConcatAA));
binCenters = edges(1:end-1) + diff(edges)/2;

figure('position',[346,123,793,255]);
subplot(121)
plot(binCenters,normCountsAA,'-o')
xlabel('Distance')
ylabel('Frequency')
title(sprintf('Median = %.1f deg',nanmedian(fieldDistConcatAA)))
box on

ylim([0 0.55])
xlim([-1 60])

[counts,edges] = histcounts(fieldDistConcatAB,binRange);
normCountsAB = counts / sum(~isnan(fieldDistConcatAB));
binCenters = edges(1:end-1) + diff(edges)/2;

subplot(122)
plot(binCenters,normCountsAB,'-o')
xlabel('Distance')
ylabel('Frequency')
title(sprintf('Median = %.1f deg',nanmedian(fieldDistConcatAB)))
box on

ylim([0 0.55])
xlim([-1 60])