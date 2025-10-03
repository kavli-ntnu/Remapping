
load expList_HP_29730
load sData_HP_29730
load ccBL_29730
sessionNums = expList_HP{1,5};

%% extract values
sessionNums = 1:size(sessionNums,2);
maps = cell(sData{1,1}.results.N,length(sessionNums));
siStore = nan(sData{1,1}.results.N,length(sessionNums));
depth = nan(sData{1,1}.results.N,length(sessionNums));
meanRate = nan(sData{1,1}.results.N,length(sessionNums));
peakRate = nan(sData{1,1}.results.N,length(sessionNums));

for iSession = sessionNums
    maps(:,iSession) = sData{1,iSession}.results.maps;
    siStore(:,iSession) = sData{1,iSession}.results.si;
    depth(:,iSession) = sData{1,iSession}.results.dpth;
    meanRate(:,iSession) = sData{1,iSession}.results.meanRate;
    peakRate(:,iSession) = sData{1,iSession}.results.peakRate;
end
    
%%
sessionComp = expList_HP{1,8};

ccStore = nan(sData{1,1}.results.N,1);
for iRow = 1:sData{1,1}.results.N
    if ~isempty(maps{iRow,sessionComp(1,1)}) && ~isempty(maps{iRow,sessionComp(1,2)})
        cc = analyses.spatialCrossCorrelation(maps{iRow,sessionComp(1,1)},maps{iRow,sessionComp(1,2)});
        ccStore(iRow,1) = cc;
    end
end

%% find place fields
spatialThresh = 0.5;
siThresh = 0.67;

fieldThresh = 0.3;
binWidth = 5;
minBins = 8;
minPeak = 1;
threshRate = 7;

numCells = size(meanRate,1);
fieldStore = nan(size(meanRate,1),2);

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
    depthMatch = repmat(1,[size(siStore,1) 1]);  
    placeIdx(:,iSession) = siStore(:,iSession) > siThresh & meanRate(:,iSession) < threshRate & meanRate(:,iSession) > 0.1 & fieldStore(:,1) > 0 & ccBaselineStability(:,iSession) > spatialThresh & depthMatch > 0;
end

placeIdx = sum(placeIdx,2);
numPlaceCells = sum(placeIdx > 0)

cellNums = (1:numCells)';
placeCells = cellNums(logical(placeIdx));

%%
selectedCells = [334 341 372 349];
figure;
cnt = 1;
for iCell = 1:size(selectedCells,2)
    for iSesh = sessionNums
        subplot(4,4,cnt)
        colorMapBRK(maps{selectedCells(iCell),iSesh},'ydir','normal');
        title(sprintf('Cell %d',selectedCells(iCell)))
        
        pk = round(peakRate(selectedCells(iCell),iSesh));
        text(0,-3,sprintf('%d Hz',pk))
        cnt = cnt + 1;
    end
end
