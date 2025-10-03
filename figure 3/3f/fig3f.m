load fig3de_sessionPairs
load expList_HP_3f

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

figure;
scatter(minNormDist,meanStore,55,idx,'filled')
close(gcf)

groupIDs = [idx(1) idx(16) idx(7)];

xVals = minNormDist;
yVals = meanStore;
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);

[r,p] = corr(xVals,yVals)

%%
pvStore = nan(size(expList_HP,1),2);
errStorePV = nan(size(expList_HP,1),2);
for iExp = 1:size(expList_HP,1)
    if ~isempty(expList_HP{iExp,9})
        load(sprintf('%s',expList_HP{iExp,8}))
        expCol = expList_HP{iExp,9};
        pvStore(iExp,1) = nanmean(simplePVstoreRaw(:,expCol,2));
        pvStore(iExp,2) = nanmedian(simplePVstoreRaw(:,expCol,2));          
        errStorePV(iExp,1) = nanstd(simplePVstoreRaw(:,expCol,2)) ./ sqrt(sum(~isnan(simplePVstoreRaw(:,expCol,2))));
    end
end

%%
load colorMat

figure;
for iGrp = 1:size(groupIDs,2)
    inds = find(idx == groupIDs(iGrp));
    
    for iExp = 1:3
        expInds = inds(expIdx(inds) == iExp);
        
        if ~isempty(expInds) && iExp == 1
            for iPt = 1:size(expInds,1)
                scatter(minNormDist(expInds(iPt),1),pvStore(expInds(iPt),2),55,colorMat(iGrp,:),'o')
                hold on
            end
        elseif ~isempty(expInds) && iExp == 2
            for iPt = 1:size(expInds,1)
                scatter(minNormDist(expInds(iPt),1),pvStore(expInds(iPt),2),55,colorMat(iGrp,:),'s')
                hold on
            end 
        elseif ~isempty(expInds) && iExp == 3
            for iPt = 1:size(expInds,1)
                scatter(minNormDist(expInds(iPt),1),pvStore(expInds(iPt),2),55,colorMat(iGrp,:),'^')
                hold on
            end
        end
    end
end

xlim([-0.02 0.72])
ylim([-0.2 0.8])
set(gca,'ytick',-0.2:0.2:1)
set(gca,'xtick',0:0.1:0.7)
xlabel('Minimum distance')
ylabel('PV correlation')
box on

xVals = minNormDist;
yVals = pvStore(:,2);
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);

[r,p] = corr(xVals,yVals)

%%
nBins = 100;

for iGrp = 1:3
    figure;
    inds = find(idx == groupIDs(iGrp));
 
    x1 = nan(size(inds,1),nBins);
    y1 = nan(size(inds,1),nBins);
    counts1 = nan(size(inds,1),nBins);
    
    for i = 1:size(inds,1)
        load(sprintf('%s',expList_HP{inds(i),8}))
        expCol = expList_HP{inds(i),9};
        temp = simplePVstoreRaw(:,expCol,2);
        
        h = histfit(temp,nBins,'kernel');
        x1(i,:) = h(2).XData;
        y1(i,:) = h(2).YData;
        counts1(i,:) = histcounts(temp,nBins);
    end
       
    xStore{1,iGrp} = nanmean(x1);
    yStore{1,iGrp} = nanmean(y1);
    countStore{1,iGrp} = nanmedian(counts1);
    
end
% close(gcf)

%%
stabMat = [];

x1 = nan(size(expList_HP,1),nBins); 
y1 = nan(size(expList_HP,1),nBins); 
counts1 = nan(size(expList_HP,1),nBins); 

for i = 1:size(expList_HP,1)
    if ~isempty(expList_HP{i,10})
        load(sprintf('%s',expList_HP{i,8}))
        expCol = expList_HP{i,10};
        temp = simplePVstoreRaw(:,expCol,2);
        stabMat = [stabMat; temp];
    
        h = histfit(temp,nBins,'kernel');
        x1(i,:) = h(2).XData;
        y1(i,:) = h(2).YData;
        counts1(i,:) = histcounts(temp,nBins);
    end
end

xStore{1,4} = nanmean(x1);
yStore{1,4} = nanmean(y1);
countStore{1,4} = nanmedian(counts1);
grpPVstore{1,4} = stabMat;

%%
figure;
for iGrp = 1:3
    currentGrp = groupIDs(iGrp);
    
    xTemp = xStore{1,iGrp};
    yTemp = yStore{1,iGrp};
    cTemp = countStore{1,iGrp};
    
    plot(xTemp,yTemp/sum(cTemp),'-')
    hold on
end

xTemp = xStore{1,4};
yTemp = yStore{1,4};
cTemp = countStore{1,4};

plot(xTemp,yTemp/sum(cTemp),'-')
hold on

ylim([0 0.03])
xlim([-0.21 1.01])

