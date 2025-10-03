load gridRateChanges_ed9d

figure('position',[-1806,192,1775,420]);

subplot(1,4,1)
allMods = [gridRateChanges{1,1}; gridRateChanges{1,2}; gridRateChanges{1,3}];
allStab = [gridRateChanges{2,1}; gridRateChanges{2,2}; gridRateChanges{2,3}];

cdfplot(allMods)
hold on
cdfplot(allStab)
axis square
grid off
title('All modules')
xlabel('Grid field rate change')
ylabel('Proportion')

subplot(1,4,2)
cdfplot(gridRateChanges{1,1})
hold on
cdfplot(gridRateChanges{2,1})
axis square
grid off
title('M1')
xlabel('Grid field rate change')
ylabel('Proportion')

subplot(1,4,3)
cdfplot(gridRateChanges{1,2})
hold on
cdfplot(gridRateChanges{2,2})
axis square
grid off
title('M2')
xlabel('Grid field rate change')
ylabel('Proportion')

subplot(1,4,4)
cdfplot(gridRateChanges{1,3})
hold on
cdfplot(gridRateChanges{2,3})
axis square
grid off
title('M3')
xlabel('Grid field rate change')
ylabel('Proportion')
