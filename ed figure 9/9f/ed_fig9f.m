load fig3de_sessionPairs
load expList_HP_3de

maxNormDist = nanmax(normDist,[],2);
meanNormDist = nanmean(normDist,2);

%%
meanStore = nan(size(expStore,1),1);
errStore = nan(size(expStore,1),1);

for iRow = 1:size(expStore,1)
    if ~isempty(expStore{iRow,5})
        meanStore(iRow,1) = nanmean(expStore{iRow,5});
        errStore(iRow,1) = nanstd(expStore{iRow,5}) ./ sqrt(sum(~isnan(expStore{iRow,5})));
    end
end

%%
load colorVec

figure('position',[-1499,285,1181,305]); 
subplot(131)
scatter(meanNormDist,meanStore,40,colorVec,'linewidth',1.5)

xlabel('Mean distance')
ylabel('Spatial correlation')
xlim([-0.02 0.82])
ylim([-0.02 0.8])
box on
set(gca,'xtick',0:0.2:1,'ytick',0:0.2:0.8)

%%
subplot(132)
scatter(maxNormDist,meanStore,40,colorVec,'linewidth',1.5)
hold on
xlabel('Max distance')
ylabel('Spatial correlation')
xlim([-0.02 1.02])
ylim([-0.02 0.8])
box on
set(gca,'xtick',0:0.2:1,'ytick',0:0.2:0.8)

%%
subplot(133)
scatter(normDist(:,1),meanStore,40,colorVec,'linewidth',1.5)
xlabel('M1-M2 distance')
ylabel('Spatial correlation')
xlim([-0.02 0.92])
ylim([-0.02 0.8])
box on
set(gca,'xtick',0:0.2:1,'ytick',0:0.2:0.8)

xVals = normDist(:,1);
yVals = meanStore;
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);

[r,p] = corr(xVals,yVals)
