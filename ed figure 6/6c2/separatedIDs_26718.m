
load expList_26718
load rotCCstore_26718
numMods = size(rotCCstore,4);
load pkStore_26718
load sData_26718
load gridIDstore_26718
load spacingStore_26718

cntCol = 1;
load colorMap
maxNumSelected = 25;
separatedIDstore = cell(1);

%%
for iExp = 1
    
    sessionNums = expList_MEC{iExp,5};

    %% set box size and nDeg
    boxSize = expList_MEC{iExp,7}(1);
    nBins = 40;
    nDeg = 3;

    %% extract values
    maps = cell(sData{1,2}.results.N,length(sessionNums));
    siStore = nan(sData{1,2}.results.N,length(sessionNums));
    depth = nan(sData{1,2}.results.N,length(sessionNums));
    meanRate = nan(sData{1,2}.results.N,length(sessionNums));
    peakRate = nan(sData{1,2}.results.N,length(sessionNums));
    grids = nan(sData{1,2}.results.N,3,length(sessionNums));
    
    for iSession = sessionNums
        
        maps(:,iSession) = sData{1,iSession}.results.maps;
        siStore(:,iSession) = sData{1,iSession}.results.si;
        depth(:,iSession) = sData{1,iSession}.results.dpth;
        meanRate(:,iSession) = sData{1,iSession}.results.meanRate;
        peakRate(:,iSession) = sData{1,iSession}.results.peakRate;
        grids(:,1:3,iSession) = sData{1,iSession}.results.grids;
    end

    %%
    sNum = sessionNums(1);
    
    currentSpacing = cell(1);
    currentDepth = cell(1);
    currentScore = cell(1);
    
    for iMod = 1:numMods
        currentMod = gridIDstore{1,iMod};
        
        currentSpacing{iExp,iMod} = spacingStore{iExp,iMod};
        currentDepth{iExp,iMod} = depth(currentMod,sNum);
        currentScore{iExp,iMod} = grids(currentMod,1,sNum);
    end

    %%
    depthMeans = nan(numMods,1);
    depthRange = nan(numMods,2);
    spacingMeans = nan(numMods,1);
    
    for iMod = 1:numMods
        depthMeans(iMod,1) = nanmean(currentDepth{iExp,iMod});
        depthRange(iMod,1) = nanmin(currentDepth{iExp,iMod});
        depthRange(iMod,2) = nanmax(currentDepth{iExp,iMod});
        
        spacingMeans(iMod,1) = nanmean(currentSpacing{iExp,iMod});
    end

    %%
    figure;
    bw = 4;
    
    for iMod = 1:numMods
        [f1,x1] = ksdensity(currentSpacing{iExp,iMod},'bandwidth',bw);
        [pks,loc,width,prom] = findpeaks(f1,x1);
        findpeaks(f1,x1,'minPeakHeight',1)
        hold on
        
        [m1,m2] = nanmax(pks);
        spacingLoc(iMod,1) = loc(m2);
        xlim([35 175])
    end
    xlabel('Spacing')
    
    %%
    for iMod = 1:numMods
        [f1,x1] = ksdensity(currentDepth{iExp,iMod});
        [pks,loc,width,prom] = findpeaks(f1,x1);
        [m1,m2] = nanmax(pks);
        depthLoc(iMod,1) = loc(m2);
    end
    
    figure;
    for iMod = 1:numMods
        scatter(currentDepth{iExp,iMod},currentSpacing{iExp,iMod},25,colorMap(iMod,:),'o','linewidth',1.5)
        hold on
    end
    plot(depthLoc,spacingLoc,'k*')

    %%
    for iMod = 1:numMods
        currentMod = gridIDstore{1,iMod};
        
        numSelected = maxNumSelected;
        if size(currentMod,1) < 25
            numSelected = size(currentMod,1);
        end
        
        d = currentDepth{iExp,iMod};
        s = currentSpacing{iExp,iMod};
        
        rescaleD = minions.rescaleData([d; depthLoc(iMod,1)],0,1);
        rescaleS = minions.rescaleData([s; spacingLoc(iMod,1)],0,1);
        
        cellDist = nan(size(currentMod,1),1);
        for iCell = 1:size(currentMod,1)
            cellDist(iCell,1) = pdist([rescaleD(iCell) rescaleS(iCell); rescaleD(end,1) rescaleS(end,1)]);
        end
        
        [s1,s2] = sort(cellDist);
        inds = s2(1:numSelected,1);
        
        separatedIDstore{iExp,iMod} = currentMod(inds);
        
        if PLOT == true
            plot(d(inds),s(inds),'ko')
            hold on
        end
    end
end

