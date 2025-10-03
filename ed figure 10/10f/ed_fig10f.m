load fig3de_sessionPairs
load expList_HP_3de
groupIDs = [1 3 2];

load colorMat

%%
load expStorePC_ed10f

figure;
selectedExp = [33 14 27 17 19 30];
cnt = 1;
for iExp = 1:size(selectedExp,2)
    x = expStorePC{selectedExp(iExp),4}; % pairDistA
    y = expStorePC{selectedExp(iExp),5}; % pairDistB
    
    pts = linspace(0, 40*sqrt(2), 80);
    N = histcounts2(y(:), x(:), pts, pts);
    
    subplot(1,size(selectedExp,2),cnt)
    imagesc(pts, pts, N);
    axis equal;
    set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
    cnt = cnt + 1;
    title(sprintf('%.2f',expStorePC{selectedExp(iExp),6}))
end

%%
load colorVec

rs = nan(size(expStorePC,1),1);
for iExp = 1:size(expStorePC,1)
    if ~isempty(expStorePC{iExp,6})
        rs(iExp,1) = expStorePC{iExp,6};
    end
end

figure;
scatter(minNormDist,rs,40,colorVec,'linewidth',1.5)
xlim([-0.02 0.72])
set(gca,'xtick',0:0.1:1)
xlabel('Minimum distance')
box on

xVals = minNormDist;
yVals = rs;
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);

[r,p] = corr(xVals,yVals)

%%
shiftStore = cell(1);
for i = 1:3
    inds = find(idx == groupIDs(i));
    
    shiftTemp = [];
    for j = 1:size(inds,1)
        temp = expStorePC{inds(j),3};
        shiftTemp = [shiftTemp; nanmedian(temp)];
    end        
    shiftStore{1,i} = shiftTemp;
end

%%
shiftCol = nan(size(expStorePC,1),1);
shiftErr = nan(size(expStorePC,1),1);

for iExp = 1:size(expStorePC,1)
    if ~isempty(expStorePC{iExp,3})
        shiftCol(iExp,1) = nanmedian(expStorePC{iExp,3});
        shiftErr(iExp,1) = nanstd(expStorePC{iExp,3}) ./ sqrt(sum(~isnan(expStorePC{iExp,3})));
    end
end
shiftCol(23,1) = nan;
shiftErr(23,1) = nan;

figure;
for iPt = 1:size(expStorePC,1)
    errorbar(minNormDist(iPt,1),shiftCol(iPt,1),shiftErr(iPt,1),'capsize',0,'color',[0.8 0.8 0.8],'linewidth',1.2)
    hold on
end
scatter(minNormDist,shiftCol,40,colorVec,'linewidth',1.5)
hold on

xlim([-0.02 0.72])
box on
ylim([0 32])
set(gca,'ytick',0:10:50)
ylabel('Place field shift (cm)')

xVals = minNormDist;
yVals = shiftCol;
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);

[r,p] = corr(xVals,yVals)

