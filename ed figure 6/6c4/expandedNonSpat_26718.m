
load expList_26718
load rotCCstore_26718
numMods = size(rotCCstore,4);
load pkStore_26718
load sData_26718
load gridIDstore_26718
load spacingStore_26718
load siCutoffs_26718

load colorMap
PLOT = true;
numSD = 2; % how many standard deviations

prctByRun = [0.1 0.2 0.3];
extraIDstore = cell(1);

for iExp = 1 %expNums
    
    sessionNums = expList_MEC{iExp,5};
    
    %% set box size and nDeg
    boxSize = expList_MEC{iExp,7}(1);
    nBins = 40;
    nDeg = 3;
    
    %% extract values
    maps = cell(sData{1,1}.results.N,length(sessionNums));
    siStore = nan(sData{1,1}.results.N,length(sessionNums));
    depth = nan(sData{1,1}.results.N,length(sessionNums));
    meanRate = nan(sData{1,1}.results.N,length(sessionNums));
    peakRate = nan(sData{1,1}.results.N,length(sessionNums));
    grids = nan(sData{1,1}.results.N,3,length(sessionNums));
    
    for iSession = sessionNums
        maps(:,iSession) = sData{1,iSession}.results.maps;
        siStore(:,iSession) = sData{1,iSession}.results.si;
        depth(:,iSession) = sData{1,iSession}.results.dpth;
        meanRate(:,iSession) = sData{1,iSession}.results.meanRate;
        peakRate(:,iSession) = sData{1,iSession}.results.peakRate;
        grids(:,1:3,iSession) = sData{1,iSession}.results.grids;
    end
    
    %%
    currentSpacing = cell(1);
    currentDepth = cell(1);
    currentScore = cell(1);
    currentSI = cell(1);
    
    for iMod = 1:numMods
        currentMod = gridIDstore{1,iMod};
        
        currentSpacing{iExp,iMod} = spacingStore{iExp,iMod};
        currentDepth{iExp,iMod} = depth(currentMod,1);
        currentScore{iExp,iMod} = grids(currentMod,1,1);
        currentSI{iExp,iMod} = siStore(currentMod,1);
    end
    
    %%
    depthMeans = nan(numMods,1);
    depthRange = nan(numMods,2);
    depthStd = nan(numMods,1);
    
    spacingMeans = nan(numMods,1);
    spacingRange = nan(numMods,2);
    spacingStd = nan(numMods,1);
    
    depthCutoffs = nan(numMods,2);
    spacingCutoffs = nan(numMods,2);
    
    for iMod = 1:numMods
        depthMeans(iMod,1) = nanmean(currentDepth{iExp,iMod});
        depthRange(iMod,1) = nanmin(currentDepth{iExp,iMod});
        depthRange(iMod,2) = nanmax(currentDepth{iExp,iMod});
        depthStd(iMod,1) = nanstd(currentDepth{iExp,iMod});
        
        spacingMeans(iMod,1) = nanmean(currentSpacing{iExp,iMod});
        spacingRange(iMod,1) = nanmin(currentSpacing{iExp,iMod});
        spacingRange(iMod,2) = nanmax(currentSpacing{iExp,iMod});
        spacingStd(iMod,1) = nanstd(currentSpacing{iExp,iMod});
    end
    
    depthCutoffs(:,1) = depthMeans - depthStd*numSD;
    depthCutoffs(:,2) = depthMeans + depthStd*numSD;
    
    spacingCutoffs(:,1) = spacingMeans - spacingStd*numSD;
    spacingCutoffs(:,2) = spacingMeans + spacingStd*numSD;
    
    %%
    depthLoc = nan(numMods,1);
    for iMod = 1:numMods
        [f1,x1] = ksdensity(currentDepth{iExp,iMod});
        [pks,loc,width,prom] = findpeaks(f1,x1);
        [m1,m2] = nanmax(pks);
        depthLoc(iMod,1) = loc(m2);
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
    rowCheck = sum(siStore > siCutoffs,2);
    selectedSI = logical(rowCheck);
    selectedSI = ~selectedSI;
    
    %%
    figure;
    currentAddition = cell(1);
    
    for iMod = 1:numMods
        currentMod = gridIDstore{1,iMod};
        
        d = depth(:,1);
        s = spacingStoreAll(:,2);
        
        rescaleD = minions.rescaleData([d; depthLoc(iMod,1)],0,1);
        rescaleS = minions.rescaleData([s; spacingLoc(iMod,1)],0,1);
        
        cellDist = nan(size(depth,1),1);
        for iCell = 1:size(depth,1)
            cellDist(iCell,1) = pdist([rescaleD(iCell) rescaleS(iCell); rescaleD(end,1) rescaleS(end,1)]);
        end
        cellDist(currentMod) = nan;
        
        intInds = logical(sum(meanRate > 10,2));
        cellDist(intInds) = nan;
        
        indsToKeep = ~logical(sum(siStore > siCutoffs,2));
        cellDist(~indsToKeep) = nan;
        
        [s1,s2] = sort(cellDist);
        
        numRuns = size(prctByRun,2);
        
        for iRun = 1:numRuns
            prct = prctByRun(iRun);
            numSelected = round(size(currentMod,1)*prct);
            inds = s2(1:numSelected,1);
            
            currentAddition{iRun,iMod} = [currentMod; inds];
            
            plot(d(inds),s(inds),'+','color',colorMap(iMod,:))
            hold on
        end
        
        % original module, box surrounds 2 std dev
        x = depthCutoffs(iMod,1);
        y = spacingCutoffs(iMod,1);
        w = depthCutoffs(iMod,2) - depthCutoffs(iMod,1);
        h = spacingCutoffs(iMod,2) - spacingCutoffs(iMod,1);
        rectangle('position',[x,y,w,h])
        
        % original module, centers
        plot(depthLoc(iMod),spacingLoc(iMod),'k*')
        
        
        for iMod = 1:numMods
            x = depthCutoffs(iMod,1);
            y = spacingCutoffs(iMod,1);
            w = depthCutoffs(iMod,2) - depthCutoffs(iMod,1);
            h = spacingCutoffs(iMod,2) - spacingCutoffs(iMod,1);
            
            rectangle('position',[x,y,w,h])
        end
    end
    extraIDstore{iExp,1} = currentAddition;
end
