load fig3de_sessionPairs
load expList_HP_3de

%%
meanStore = nan(size(expStore,1),1);
errStore = nan(size(expStore,1),1);

for iRow = 1:size(expStore,1)
    if ~isempty(expStore{iRow,5})
        meanStore(iRow,1) = nanmean(expStore{iRow,5});
        errStore(iRow,1) = nanstd(expStore{iRow,5}) ./ sqrt(sum(~isnan(expStore{iRow,5})));
    end
end

meanStoreStab = nan(size(expStore,1),1);
errStoreStab = nan(size(expStore,1),1);

for iRow = 1:size(expStoreStab,1)
    if ~isempty(expStoreStab{iRow,5})
        meanStoreStab(iRow,1) = nanmean(expStoreStab{iRow,5});
        errStoreStab(iRow,1) = nanstd(expStoreStab{iRow,5}) ./ sqrt(sum(~isnan(expStoreStab{iRow,5})));
    end
end

%%
idx = kmeans([minNormDist meanStore],3,'replicates',1000);
groupIDs = [idx(1) idx(16) idx(7)];

figure;
scatter(minNormDist,meanStore,55,idx,'filled')
close(gcf)

%%
xVals = minNormDist;
yVals = meanStore;
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);

[r,p] = corr(xVals,yVals)

%%
load colorMat

figure;
for iGrp = 1:size(groupIDs,2)
    inds = find(idx == groupIDs(iGrp));
    
    for iExp = 1:3
        expInds = inds(expIdx(inds) == iExp);
        
        if ~isempty(expInds) && iExp == 1
            for iPt = 1:size(expInds,1)
                scatter(minNormDist(expInds(iPt),1),meanStore(expInds(iPt),1),55,colorMat(iGrp,:))
                hold on
                errorbar(minNormDist(expInds(iPt),1),meanStore(expInds(iPt),1),errStore(expInds(iPt),1),'capsize',0,'color',[0.8 0.8 0.8],'linewidth',1.2)
            end
        elseif ~isempty(expInds) && iExp == 2
            for iPt = 1:size(expInds,1)
                scatter(minNormDist(expInds(iPt),1),meanStore(expInds(iPt),1),55,colorMat(iGrp,:),'s')
                hold on
                errorbar(minNormDist(expInds(iPt),1),meanStore(expInds(iPt),1),errStore(expInds(iPt),1),'-','capsize',0,'color',[0.8 0.8 0.8],'linewidth',1.2)
            end 
        elseif ~isempty(expInds) && iExp == 3
            for iPt = 1:size(expInds,1)
                scatter(minNormDist(expInds(iPt),1),meanStore(expInds(iPt),1),55,colorMat(iGrp,:),'^')
                errorbar(minNormDist(expInds(iPt),1),meanStore(expInds(iPt),1),errStore(expInds(iPt),1),'-','capsize',0,'color',[0.8 0.8 0.8],'linewidth',1.2)
                hold on
            end
        end
    end
end

xlim([-0.02 0.72])
ylim([-0.05 0.72])
set(gca,'ytick',-0.2:0.1:1)
set(gca,'xtick',0:0.1:0.7)
xlabel('Minimum distance')
ylabel('Spatial correlation')
box on

%%
idx1 = find(idx == 1);
grp1 = [];
for i = 1:size(idx1,1)
    grp1 = [grp1; expStore{idx1(i),5}];
end
   
figure;
h1 = histfit(grp1,100,'kernel');
x1 = h1(2).XData;
y1 = h1(2).YData;
counts1 = histcounts(grp1,100);

close(gcf)

idx2 = find(idx == 2);
grp2 = [];
for i = 1:size(idx2,1)
    grp2 = [grp2; expStore{idx2(i),5}];
end

figure;
h2 = histfit(grp2,100,'kernel');
x2 = h2(2).XData;
y2 = h2(2).YData;
counts2 = histcounts(grp2,100);

close(gcf)

idx3 = find(idx == 3);
grp3 = [];
for i = 1:size(idx3,1)
    grp3 = [grp3; expStore{idx3(i),5}];
end

figure;
h3 = histfit(grp3,100,'kernel');
x3 = h3(2).XData;
y3 = h3(2).YData;
counts3 = histcounts(grp3,100);

close(gcf)

%%
stabGrp = [];
for iRow = 1:size(expStoreStab,1)
    if ~isempty(expStoreStab{iRow,5})
        temp = expStoreStab{iRow,5};
        stabGrp = [stabGrp; temp];
    end
end
   
figure;
h4 = histfit(stabGrp,100,'kernel');
x4 = h4(2).XData;
y4 = h4(2).YData;
counts4 = histcounts(stabGrp,100);

close(gcf)

%% plot them together!
figure;
plot(x1,y1/sum(counts1),'-')
hold on
plot(x2,y2/sum(counts2),'-')
plot(x3,y3/sum(counts3),'-')
plot(x4,y4/sum(counts4),'-')

xlabel('Spatial correlation')
ylabel('Frequency')
ylim([0 0.03])
xlim([-0.52 1.02])


