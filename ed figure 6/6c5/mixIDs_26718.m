
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
prctByRun = [0.1 0.3 0.5];

mixIDstore = cell(1);

for iExp = 1 %expNums
    
    if numMods == 3
        compStore = [2 nan; 1 3; 2 nan];
    else
        compStore = [2 nan; 1 nan];
    end

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
    
    for iMod = 1:numMods
        currentMod = gridIDstore{1,iMod};
        
        currentSpacing{iExp,iMod} = spacingStore{iExp,iMod};
        currentDepth{iExp,iMod} = depth(currentMod,1);
        currentScore{iExp,iMod} = grids(currentMod,1,1);
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
    currentMix = cell(1);
    
    numPlts = nansum(nansum(~isnan(compStore)));
    colCnt = 1;
    
    if PLOT == true
        figure('position',[-1780,130,1758,320]);
    end
    
    for iMod = 1:numMods
        
        currentMod = gridIDstore{1,iMod};
        kpInds = ~isnan(compStore(iMod,:));
        modsToComp = compStore(iMod,kpInds);
        
        for iComp = 1:size(modsToComp,2)
            
            if PLOT == true
                subplot(1,numPlts,colCnt)
                for i = 1:numMods
                    scatter(currentDepth{iExp,i},currentSpacing{iExp,i},25,colorMap(i,:),'o','linewidth',1.5)
                    hold on
                    plot(depthLoc,spacingLoc,'k*')
                end
            end
            
            compMod = gridIDstore{1,modsToComp(iComp)};
            
            d = currentDepth{iExp,modsToComp(iComp)};
            s = currentSpacing{iExp,modsToComp(iComp)};
            
            rescaleD = minions.rescaleData([d; depthLoc(iMod,1)],0,1);
            rescaleS = minions.rescaleData([s; spacingLoc(iMod,1)],0,1);
            
            cellDist = nan(size(d,1),1);
            for iCell = 1:size(d,1)
                cellDist(iCell,1) = pdist([rescaleD(iCell) rescaleS(iCell); rescaleD(end,1) rescaleS(end,1)]);
            end
            
            [s1,s2] = sort(cellDist);
            
            numRuns = size(prctByRun,2);
            
            for iRun = 1:numRuns
                prct = prctByRun(iRun);
                numSelected = round(size(currentMod,1)*prct);
                
                if numSelected > size(compMod,1)
                    numSelected = size(compMod,1);
                end
                inds = s2(1:numSelected,1);
                
                currentMix{iRun,colCnt} = [currentMod; compMod(inds)];
                
                if PLOT == true
                    plot(d(inds),s(inds),'ko')
                    hold on
                    title(sprintf('M%d + M%d',iMod,modsToComp(iComp)))
                    xlabel('Distance to tip of probe')
                    ylabel('Spacing')
                end
            end
            colCnt = colCnt + 1;
        end
    end
    mixIDstore{iExp,1} = currentMix;
end

