load expList_26718
load rotCCstore_26718
numMods = size(rotCCstore,4);
load pkStore_26718
load singleCells_26718

cntCol = 1;

%%
fieldDistStore = cell(1);
otherDistStore = cell(1);

for iExp = 1
    for iMod = 1:numMods
        ccStoreAll = singleCells{1,iMod};
        units = 1:size(ccStoreAll,4);
        
        modNums = 1:numMods;
        modNums(iMod) = nan;
        modNums = modNums(~isnan(modNums));
        
        fieldDist = nan(size(units,2),1);
        otherDist = nan(size(units,2),size(modNums,2));
        for iCell = 1:size(units,2)
            
            if ~isempty(ccStoreAll(:,:,mindStore(iMod),iCell))
                cc = ccStoreAll(:,:,mindStore(iMod),iCell);
                
                cc(isnan(cc)) = 0;
                rm = imregionalmax(cc);
                rp = regionprops(rm);
                numPeaks = size(rp,1);
                
                vx = []; vy = [];
                for i = 1:length(rp)
                    vx = [vx; rp(i).Centroid(1)];
                    vy = [vy; rp(i).Centroid(2)];
                end
                
                [idx,d] = knnsearch([vx vy],pkStore(iMod,:)); % find the points closest to the pk for each cell
                fieldDist(iCell,1) = d;
                
                cnt = 1;
                for i = modNums
                    otherDist(iCell,cnt) = pdist([vx(idx) vy(idx); pkStore(i,1) pkStore(i,2)]);
                    cnt = cnt + 1;
                end
            end
        end
        
        nanmedian(fieldDist)
        nanmedian(otherDist)
        
        fieldDistStore{1,cntCol} = fieldDist;
        otherDistStore{1,cntCol} = otherDist;
        cntCol = cntCol + 1;
    end
end

