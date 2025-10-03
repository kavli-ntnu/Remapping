load expList_HP_29730
load sData_HP_29730
load placeCells_29730
numCells = size(placeCells,1);

iExp = 1;
sessionNums = expList_HP{iExp,5};
numSessions = size(sessionNums,2);
sessionComp = [1 2];

depthFilter = [];

sesh1 = sessionComp(1);
sesh2 = sessionComp(2);

%% extract values
maps = cell(sData{1,1}.results.N,length(sessionNums));
for iSession = sessionNums
    maps(:,iSession) = sData{1,iSession}.results.maps;
end

%%
dpStore = [];
for iCell = 1:size(placeCells,1)
    mapA = maps{placeCells(iCell),sesh1};
    mapB = maps{placeCells(iCell),sesh2};
    
    [r1,c1] = find(mapA == nanmax(nanmax(mapA)));
    [r2,c2] = find(mapB == nanmax(nanmax(mapB)));
    
    if ~isempty(c1) & ~isempty(r1) & ~isempty(c2) & ~isempty(r2)
        p1 = [c1(1) r1(1)]; 
        p2 = [c2(1) r2(1)];
        dp = p2-p1; % Difference
        dpStore(iCell,:) = dp;
    end
end

%%
entries = unique(dpStore,'rows');
uniqueEntries = minions.removeNans(entries,'rows','any');

numEntries = [];
for iRow = 1:size(uniqueEntries,1)
    temp = uniqueEntries(iRow,:);
    numEntries(iRow,1) = sum(dpStore(:,1) == temp(1) & dpStore(:,2) == temp(2));
end
maxNum = nanmax(numEntries);

expStore{iExp,1} = expList_HP{iExp,3};
expStore{iExp,2} = sessionComp;
expStore{iExp,3} = dpStore;
expStore{iExp,4} = uniqueEntries;
expStore{iExp,5} = numEntries;

%% store peaks
pkStore = nan(numCells,4);
distStore = nan(numCells,1);

for iCell = 1:size(placeCells,1)
    
    mapA = maps{placeCells(iCell),sesh1};
    mapB = maps{placeCells(iCell),sesh2};
    
    [r1,c1] = find(mapA == nanmax(nanmax(mapA)));
    [r2,c2] = find(mapB == nanmax(nanmax(mapB)));
    
    if ~isempty(r1) && ~isempty(c1) && ~isempty(r2) && ~isempty(c2)
        pkStore(iCell,1:4) = [c1(1) r1(1) c2(1) r2(1)];
        distStore(iCell,1) = pdist([c1(1) r1(1); c2(1) r2(1)]);
    end
end

expStore{iExp,6} = minions.removeNans(pkStore,'rows','all');
expStore{iExp,7} = minions.removeNans(distStore,'rows','all');

%% remapping val scatter
storeCnt = 1;
cnt = 1;

pairDistA = nan(15000,1);
pairDistB = nan(15000,1);

pairRateA = nan(15000,1);
pairRateB = nan(15000,1);

idx = ~isnan(pkStore(:,1));
rowIDs = (1:numCells)';
allCells = (rowIDs(idx))';

for iCell = allCells
    
    mapA = maps{placeCells(iCell),sesh1};
    mapB = maps{placeCells(iCell),sesh2};
    
    idx = ~isnan(pkStore(:,1));
    rowIDs = (1:numCells)';
    cellIDs = rowIDs(idx);
    
    currentA = pkStore(iCell,1:2);
    currentB = pkStore(iCell,3:4);
    
    currentRateA = nanmax(nanmax(mapA));
    currentRateB = nanmax(nanmax(mapB));
    
    cellIDs(cnt) = nan;
    idsToCheck = (cellIDs(~isnan(cellIDs)))';
    
    for j = idsToCheck
        otherA = pkStore(j,1:2);
        otherB = pkStore(j,3:4);
        
        mapOtherA = maps{placeCells(j),sesh1};
        mapOtherB = maps{placeCells(j),sesh2};
        
        otherRateA = nanmax(nanmax(mapOtherA));
        otherRateB = nanmax(nanmax(mapOtherB));
        
        tempDistA = pdist([currentA; otherA]);
        tempDistB = pdist([currentB; otherB]);
        
        tempRateA = calc.diffScore(currentRateA,otherRateA);
        tempRateB = calc.diffScore(currentRateB,otherRateB);
        
        pairDistA(storeCnt,1) = tempDistA;
        pairDistB(storeCnt,1) = tempDistB;
        
        pairRateA(storeCnt,1) = tempRateA;
        pairRateB(storeCnt,1) = tempRateB;
        
        storeCnt = storeCnt + 1;
    end
    cnt = cnt + 1;
end

n = sum(~isnan(pairDistA))

pairDistA = pairDistA(~isnan(pairDistA));
pairDistB = pairDistB(~isnan(pairDistA));

[r,p] = corr(pairDistA,pairDistB)
remappingVal = 1-r

expStore{iExp,8} = pairDistA;
expStore{iExp,9} = pairDistB;
expStore{iExp,10} = remappingVal;

%% remapping val scatter and heat map
x = pairDistA;
y = pairDistB;

pts = linspace(0, 40*sqrt(2), 80);
N = histcounts2(y(:), x(:), pts, pts);

% Plot scattered data (for comparison):
figure;
imagesc(pts, pts, N);
axis equal;
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');


