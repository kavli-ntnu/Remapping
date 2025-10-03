load expList_HP_29730
load sData_HP_29730
load ccBL_29730

iExp = 1;
sessionNums = expList_HP{iExp,5};
sessionComp = expList_HP{iExp,8};

%%
maps = cell(sData{1,1}.results.N,length(sessionNums));
siStore = nan(sData{1,1}.results.N,length(sessionNums));
depth = nan(sData{1,1}.results.N,length(sessionNums));
meanRate = nan(sData{1,1}.results.N,length(sessionNums));

for iSession = sessionNums
    maps(:,iSession) = sData{1,iSession}.results.maps;
    siStore(:,iSession) = sData{1,iSession}.results.si;
    depth(:,iSession) = sData{1,iSession}.results.dpth;
    meanRate(:,iSession) = sData{1,iSession}.results.meanRate;
end

%%
ccStore = nan(sData{1,1}.results.N,1);
for iRow = 1:sData{1,1}.results.N
    if ~isempty(maps{iRow,1}) && ~isempty(maps{iRow,4})
        cc = analyses.spatialCrossCorrelation(maps{iRow,sessionNums(1)},maps{iRow,sessionNums(4)});
        ccStore(iRow,1) = cc;
    end
end

%% find place fields 
numCells = size(meanRate,1);
fieldStore = nan(size(meanRate,1),2);

% set thresholds
spatialThresh = 0.5;
siThresh = 0.67;

fieldThresh = 0.3;
binWidth = 5;
minBins = 8;
minPeak = 1;
threshRate = 7;

for iSession = sessionNums
    for iRow = 1:numCells
        [~,fields] = analyses.placefield(maps{iRow,iSession},'threshold',fieldThresh,'binWidth',binWidth,'minBins',minBins,'minPeak',minPeak);
        if ~isempty(fields)
            fieldStore(iRow,iSession) = length(fields);
        end
    end
end

placeIdx = [];
for iSession = sessionNums  
    placeIdx(:,iSession) = siStore(:,iSession) > siThresh & meanRate(:,iSession) < threshRate & meanRate(:,iSession) > 0.1 & fieldStore(:,1) > 0 & ccBaselineStability(:,iSession) > 0.5; 
end

placeIdx = sum(placeIdx,2);
numPlaceCells = sum(placeIdx > 0)

cellNums = (1:numCells)';
placeCells = cellNums(logical(placeIdx));

%% spatial correlation
placeCellCorr = nan(size(placeCells,1),size(sessionComp,1));

for iComp = 1:size(sessionComp,1)
    for iRow = 1:size(placeCells,1)
        if ~isempty(maps{placeCells(iRow),sessionComp(iComp,1)}) && ~isempty(maps{placeCells(iRow),sessionComp(iComp,2)})
            cc = analyses.spatialCrossCorrelation(maps{placeCells(iRow),sessionComp(iComp,1)},maps{placeCells(iRow),sessionComp(iComp,2)});
            placeCellCorr(iRow,iComp) = cc;
        end
    end
end
nanmean(placeCellCorr)

%%
placeCellDepth = depth(placeCells,1);

%%
placeRates = nan(size(placeCells,1),size(sessionNums,2));
for iRow = 1:size(placeCells,1)
    for iSession = sessionNums
        placeRates(iRow,iSession) = meanRate(placeCells(iRow),iSession);
    end
end

rateDiff = nan(size(placeCells,1),size(sessionComp,1));
for iComp = 1:size(sessionComp,1)
    rateDiff(:,iComp) = abs(calc.diffScore(placeRates(:,sessionComp(iComp,1)),placeRates(:,sessionComp(iComp,2))));
end
nanmean(rateDiff)
    