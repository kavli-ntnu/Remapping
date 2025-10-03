load fig3de_sessionPairs
load expList_HP_3de
load spacingDiff_ed9

%%
meanStore = nan(size(expStore,1),1);

for iRow = 1:size(expStore,1)
    if ~isempty(expStore{iRow,5})
        meanStore(iRow,1) = nanmean(expStore{iRow,5});
    end
end

%%
load colorVec

figure;
scatter(spacingDiff,meanStore,40,colorVec,'linewidth',1.5)
xlim([0 85])
ylim([-0.02 0.8])
set(gca,'ytick',0:0.2:1)
box on
xlabel('Spacing difference (cm)')
ylabel('Spatial correlation')

xVals = meanStore;
yVals = spacingDiff;
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);

[r,p] = corr(xVals,yVals)

%%
figure;
scatter(spacingDiff,minNormDist,40,colorVec,'linewidth',1.5)
hold on
xlim([0 85])
ylim([-0.02 0.8])
set(gca,'ytick',0:0.2:1)
box on
xlabel('Spacing difference (cm)')
ylabel('Minimum distance')

xVals = spacingDiff;
yVals = minNormDist;
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);

[r,p] = corr(xVals,yVals)
