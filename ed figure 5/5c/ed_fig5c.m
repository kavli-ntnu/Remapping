%%
load cellPairDist_AA_ed5c
maxValAA = nanmax([nanmax(withinDistColAA) nanmax(betweenDistColAA)]);

binRange = 0:10:round(maxValAA);
[counts,edges] = histcounts(withinDistColAA,binRange);
normCountsWithinAA = counts / sum(~isnan(withinDistColAA));
binCenters = edges(1:end-1) + diff(edges)/2;

figure;
hold on
plot(binCenters,normCountsWithinAA,'-o','color',[0.5 0.5 0.5])

[counts,edges] = histcounts(betweenDistColAA,binRange);
normCountsBetweenAA = counts / sum(~isnan(betweenDistColAA));
binCenters = edges(1:end-1) + diff(edges)/2;

hold on
plot(binCenters,normCountsBetweenAA,'-o','color','k')
xlabel('Distance')
ylabel('Frequency')
ylim([0 1])
xlim([0 60])
set(gca,'ytick',0:0.2:1)
box on

%%
load cellPairDist_AB_ed5c

maxValAB = nanmax([nanmax(withinDistColAB) nanmax(betweenDistColAB)]);

binRange = 0:10:round(maxValAB);
[counts,edges] = histcounts(withinDistColAB,binRange);
normCountsWithinAB = counts / sum(~isnan(withinDistColAB));
binCenters = edges(1:end-1) + diff(edges)/2;

figure;
plot(binCenters,normCountsWithinAB,'-o')
hold on

[counts,edges] = histcounts(betweenDistColAB,binRange);
normCountsBetweenAB = counts / sum(~isnan(betweenDistColAB));
binCenters = edges(1:end-1) + diff(edges)/2;

plot(binCenters,normCountsBetweenAB,'k-o')

xlabel('Distance')
ylabel('Frequency')
ylim([0 0.5])
xlim([0 160])
set(gca,'ytick',0:0.1:1,'xtick',0:50:150)
box on

