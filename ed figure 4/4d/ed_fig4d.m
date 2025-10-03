load minDiffStore_ed4d

ABstore = [];
AAstore = [];
for iCol = 1:size(minDiffStore,2)
    ABstore = [ABstore; minDiffStore{1,iCol}];
    AAstore = [AAstore; minDiffStore{2,iCol}];
end

maxVal = 40;
binRange = nanmin(AAstore(:,1)):3:maxVal;

[counts,edges] = histcounts(AAstore(:,1),binRange);
normCountsAA = counts / sum(~isnan(AAstore(:,1)));
binCenters = edges(1:end-1) + diff(edges)/2;

figure('position',[346,123,793,255]);
subplot(121)
plot(binCenters,normCountsAA,'-o')
xlabel('Rotation (deg)')
ylabel('Frequency')
title(sprintf('Median = %.1f deg',nanmedian(AAstore(:,1))))
box on

[counts,edges] = histcounts(ABstore(:,1),binRange);
normCountsAB = counts / sum(~isnan(ABstore(:,1)));
binCenters = edges(1:end-1) + diff(edges)/2;

subplot(122)
plot(binCenters,normCountsAB,'-o')
xlabel('Rotation (deg)')
ylabel('Frequency')
title(sprintf('Median = %.1f deg',nanmedian(ABstore(:,1))))
box on

ylim([0 0.6])
