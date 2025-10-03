
load('Q:\chrislyk\from ben\Lisman super sim\rotationsToTest.mat')
load('Q:\chrislyk\simulationPV\rotationPV_all.mat')

%%
% module rotations
for iRow = 1:size(rotationsToTest,1)
    temp = rotationsToTest{iRow,1};
    m1(iRow,1) = temp(1);
    m2(iRow,1) = temp(2);
    m3(iRow,1) = temp(3);
end

%%
minDiffMat = nan(size(rotationsToTest,1),3);
minDiffMat(:,1) = abs(m1-m2);
minDiffMat(:,2) = abs(m1-m3);
minDiffMat(:,3) = abs(m2-m3);

minDiff = nanmin(minDiffMat,[],2);

%%
compVar = maxCorrAll;

bins = 0:30;
store = nan(1,size(bins,2));
for i = 1:size(bins,2)
    store(i) = nanmean(compVar(minDiff == bins(i)));
end

%% curve fitting toolbox
xVals = bins;
yVals = store;

createFit3(xVals,yVals)

%%
set(gcf,'position',[-718,353,552,422],'color','w')
grid off
legend off
set(gca,'xtick',0:5:30)
set(gca,'ytick',0.2:0.1:1)
xlabel('Minimum difference')
ylabel('PV correlation')
xlim([-1 31])
ylim([0.2 0.8])

% add error bars
hold on
for i = 1:size(bins,2)
    temp = compVar(minDiff == bins(i));
    
    errorbar(bins(i),nanmean(store(i)),(nanstd(temp)) ./ (sqrt(sum(~isnan(temp)))),'.','markersize',10,'markeredgecolor','k','markerfacecolor','k','linewidth',1,'color','k');
    hold on
end

