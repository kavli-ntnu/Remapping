load fig3de_sessionPairs
load expList_HP_3de

%%
figure('position',[723,214,274,314]);

for iExp = 1:size(normDist,1)  
    temp = normDist(iExp,:);
    
    plot(1:size(temp,2),temp,'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:size(temp,2),temp,'o','color',[0.5 0.5 0.5],'linewidth',1.5)   
end

xlim([0.5 3.5])
ylim([-0.02 1.02])
set(gca,'xtick',1:3,'xticklabels',{'M1-M2','M1-M3','M2-M3'})
set(gca,'ytick',0:0.2:2)
ylabel('Normalized distance')
box on

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
load minDisp_ed9c

figure;
scatter(minNormDist,minDisp,50,meanStore,'filled')
box on
xlabel('Minimum distance')
ylabel('Minimum displacement')
set(gca,'xtick',0:0.1:1,'ytick',0:0.1:0.4)
xlim([-0.02 0.72])
ylim([-0.02 0.42])
% axis square
colormap(flipud(colormap))

%%
load colorVec

figure; 
scatter(minDisp,meanStore,40,colorVec,'linewidth',1.5)
hold on

xlim([-0.02 0.42])
ylim([-0.02 0.8])
set(gca,'ytick',-0.2:0.2:1)
set(gca,'xtick',0:0.1:0.5)
xlabel('Minimum displacement')
ylabel('Spatial correlation')
box on

xVals = minDisp;
yVals = meanStore;
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);

[r,p] = corr(xVals,yVals)