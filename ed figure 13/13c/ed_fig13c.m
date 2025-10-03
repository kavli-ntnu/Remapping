load heatMaps_ed13c

compVar = maxCorrAll;
matStore = nan(4,3);
numStore = nan(4,3);
rowCnt = 1;
for i = 0:3 % numMods
    colCnt = 1;
    for j = [0 1 2 3] % numPairs
        
        inds = find(numAboveCutoff == i & pairsAboveCutoff == j);
        temp = compVar(inds);
        
        cellStore{rowCnt,colCnt} = temp;
        matStore(rowCnt,colCnt) = nanmedian(temp);
        numStore(rowCnt,colCnt) = size(temp,1);
        colCnt = colCnt + 1;
    end
    rowCnt = rowCnt + 1;
end
matStorePV = matStore;

%%
compVar = rValAll;
matStore = nan(4,3);
numStore = nan(4,3);
rowCnt = 1;
for i = 0:3 % numMods
    colCnt = 1;
    for j = [0 1 2 3] % numPairs
        
        inds = find(numAboveCutoff == i & pairsAboveCutoff == j);
        temp = compVar(inds);
        
        cellStore{rowCnt,colCnt} = temp;
        matStore(rowCnt,colCnt) = nanmedian(temp);
        numStore(rowCnt,colCnt) = size(temp,1);
        colCnt = colCnt + 1;
    end
    rowCnt = rowCnt + 1;
end
matStoreRS = matStore;

%%
compVar = validStoreAll;
matStore = nan(4,3);
numStore = nan(4,3);
rowCnt = 1;
for i = 0:3 % numMods
    colCnt = 1;
    for j = [0 1 2 3] % numPairs
        
        inds = find(numAboveCutoff == i & pairsAboveCutoff == j);
        temp = compVar(inds);
        
        cellStore{rowCnt,colCnt} = temp;
        matStore(rowCnt,colCnt) = nanmedian(temp);
        numStore(rowCnt,colCnt) = size(temp,1);
        colCnt = colCnt + 1;
    end
    rowCnt = rowCnt + 1;
end
matStorePrc = matStore;

%% heat maps
figure;
colorMapBRK(flipud(matStorePV))
colormap(jet)
caxis([0.3 0.9])
title('PV correlation')

figure;
colorMapBRK(flipud(matStoreRS))
colormap(flipud(jet))
caxis([0.4 0.95])
title('Rearrangement score')

figure;
colorMapBRK(flipud(matStorePrc))
colormap(flipud(jet))
caxis([0.2 0.55])
title('Percent turnover')
