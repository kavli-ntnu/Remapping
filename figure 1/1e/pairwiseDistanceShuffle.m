
load expList_26821
load gridIDstore_26821
numMods = size(gridIDstore,2);
load rotCCstore_AB_26821

numRuns = 1000;
PLOT = false;

finalModPtStore = cell(1);
finalDistStore = cell(1);
spacingAllExpt = nan(size(expList_MEC,1),3);

%%
for iExp = 1 %expNums
    
    exptSelected = iExp;
    
    name = expList_MEC{exptSelected,3};
    boxSize = expList_MEC{exptSelected,7}(1);
    szmap = 40;
    
    selectedMods = 1:numMods;
     
    %%
    pkStore = nan(numMods,2);
    maskStore = nan(size(rotCCstore,1),size(rotCCstore,1),numMods);
    ccCoordStore = nan(7,2,numMods);
    spacingStore = nan(1,numMods);
    orientStore = nan(numMods,3);
    
    %%
    for iMod = 1:numMods
        cc = rotCCstore(:,:,mindStore(iMod),iMod);
        %         cc = rotCCstore(:,:,iMod);
        
        cc(isnan(cc)) = 0;
        ac = xcorr2(cc);
        
        [gScore,gStats,~,COM] = analyses.gridnessScore(ac);
        origin = round(size(ac,1) / 2);
        
        if size(COM,1) < 6
            cc = rotCCstore(:,:,mindStore(iMod)-1,iMod);
            cc(isnan(cc)) = 0;
            ac = xcorr2(cc);
            [gScore,gStats,~,COM] = analyses.gridnessScore(ac);
        end
        
        if isempty(gStats.spacing)
            rm = imregionalmax(ac);
            rp = regionprops(rm);
            numPeaks = size(rp,1);
            prct = 0.2;
            
            while numPeaks < 6
                threshVal = nanmax(nanmax(cc))*prct;
                ac = xcorr2(double(im2bw(cc,threshVal)));
                
                rm = imregionalmax(ac);
                rp = regionprops(rm);
                numPeaks = size(rp,1);
                prct = prct + 0.05;
            end
            
            % vx and vy coordinates from ac
            vx = []; vy = [];
            for i = 1:length(rp)
                vx = [vx; rp(i).Centroid(1)];
                vy = [vy; rp(i).Centroid(2)];
            end
            % nearest vx and vy to center of ac
            kind = knnsearch([vx,vy],[size(cc,1),size(cc,2)],'K',7);
            acX = vx(kind);
            acY = vy(kind);
            
            spacing = [];
            for iPt = 1:size(acX,1)
                spacing(iPt,1) = pdist([origin origin; acX(iPt) acY(iPt)]);
            end
            spacing = spacing(2:end);
            hLength = nanmean(spacing)*tan(pi/6);
            spacingStore(1,iMod) = nanmean(spacing)*boxSize/szmap;
            
            COM = [acX acY];
            COM = COM(2:end,:);
            
            orientation = nan(size(COM,1),1);
            for iOrient = 1:size(COM,1)
                xChange = COM(iOrient,1) - origin;
                yChange = COM(iOrient,2) - origin;
                orientation(iOrient,1) = rad2deg(atan(yChange/xChange));
            end
            orientation = unique(orientation);
            hOrient = orientation + 30;
            orientStore(iMod,1:3) = (orientation)';
        else
            
            hLength = nanmean(gStats.spacing)*tan(pi/6);
            hOrient = gStats.orientation + 30;
            spacingStore(1,iMod) = nanmean(gStats.spacing)*boxSize/szmap;
            orientStore(iMod,1:3) = (gStats.orientation)';
           
        end
        
        coords = nan(6,2);
        for iOrient = 1:3
            temp = hOrient(iOrient);
            
            x = hLength*cos(deg2rad(temp));
            y = hLength*sin(deg2rad(temp));
            
            coords(iOrient,1) = origin + x;
            coords(iOrient,2) = origin + y;
        end
        
        for iOrient = 1:3
            temp = hOrient(iOrient)+180;
            
            x = hLength*cos(deg2rad(temp));
            y = hLength*sin(deg2rad(temp));
            
            coords(iOrient+3,1) = origin + x;
            coords(iOrient+3,2) = origin + y;
        end
        
        ccCoords = coords - 39;
        
        figure;
        subplot(1,3,1)
        colorMapBRK(ac);
        axis on
        
        hold on
        plot(COM(:,1),COM(:,2),'k.')
        plot(coords(:,1),coords(:,2),'wo')
        
        subplot(1,3,2)
        colorMapBRK(rotCCstore(:,:,mindStore(iMod),iMod));
        %     colorMapBRK(rotCCstore(:,:,iMod));
        hold on
        axis on
        xline(39)
        yline(39)
        plot(ccCoords(:,1),ccCoords(:,2),'k.','markersize',25)
        
        vx = ccCoords(:,1);
        vy = ccCoords(:,2);
        
        ccWrapped = vertcat(ccCoords,ccCoords(1,:));
        ccCoordStore(:,:,iMod) = ccWrapped;
        
        % find peak in xcorr
        [X,Y] = meshgrid(1:size(cc,1),1:size(cc,2));
        % which points X and Y are within polygon defined by vx + vy
        % this creates a mask to overlay onto cc, turns everything
        % outside polygon to 0
        idx = inpolygon(X,Y,vx,vy);
        xcmask = cc;
        xcmask(~idx) = nan;
        % find the highest value in the polygon
        [~,mindRhom] = max(xcmask,[],'all','linear');
        [y,x] = ind2sub([size(cc,1),size(cc,2)],mindRhom);
        
        subplot(1,3,3);
        set(gcf,'color','w')
        colorMapBRK(xcmask);
        axis on
        hold on
        plot(x,y,'ko','markers',15)
        
        pkStore(iMod,1) = x;
        pkStore(iMod,2) = y;
        
        maskStore(:,:,iMod) = xcmask;
        
%         close(gcf)
    end
    
    
    %%
    load colorMap
    load colorMapLight
    comp = [2 4; 4 2; 6 4; 4 6; 5 3; 3 5];
    
    numReps = 20;
    modPks = cell(1);
    
    for iMod = 1:numMods
        figure('position',[1,41,1280,683]);
        
        storeAllPeaks = cell(1);
        
        vx = ccCoordStore(:,1,iMod);
        vy = ccCoordStore(:,2,iMod);
        tileStore{1,1} = [vx vy];
        
        plot(vx,vy,'k-')
        
        changeStore = nan(size(comp,1),2);
        % find all xChange and yChange values
        for i = 1:size(comp,1)
            xChange = vx(comp(i,1))-vx(comp(i,2));
            yChange = vy(comp(i,1))-vy(comp(i,2));
            
            changeStore(i,1) = xChange;
            changeStore(i,2) = yChange;
        end
        
        hold on
        %     plot(pkStore(iMod,1),pkStore(iMod,2),'o','color',colorMap(iMod,:),'linewidth',2,'markerfacecolor',colorMapLight(iMod,:))
        plot(pkStore(iMod,1),pkStore(iMod,2),'o','color',colorMap(iMod,:),'linewidth',2)
        xDist = vx(1)-pkStore(iMod,1);
        yDist = vy(1)-pkStore(iMod,2);
        
        % tile center tile up
        [~,~,peakStore] = tileUp(vx,vy,changeStore,numReps,xDist,yDist);
        storeAllPeaks{1,1} = peakStore;
        
        axis square
        
        % tile center tile down
        vx = ccCoordStore(:,1,iMod);
        vy = ccCoordStore(:,2,iMod);
        [~,~,peakStore] = tileDown(vx,vy,changeStore,numReps,xDist,yDist);
        storeAllPeaks{2,1} = peakStore;
        
        %%
        vx = ccCoordStore(:,1,iMod);
        vy = ccCoordStore(:,2,iMod);
        
        [currentTileX,currentTileY,peakStore] = tileChange(vx,vy,3,changeStore,numReps,xDist,yDist);
        rightTile = [currentTileX currentTileY];
        
        [currentTileX,currentTileY,peakStore] = tileChange(vx,vy,4,changeStore,numReps,xDist,yDist);
        leftTile = [currentTileX currentTileY];
        
        xChangeR = rightTile(:,1) - tileStore{1}(:,1);
        xChangeL = leftTile(:,1) - tileStore{1}(:,1);
        
        yChangeR = rightTile(:,2) - tileStore{1}(:,2);
        yChangeL = leftTile(:,2) - tileStore{1}(:,2);
        
        %% RIGHT
        currentTileX = rightTile(:,1);
        currentTileY = rightTile(:,2);
        currentPk(1) = currentTileX(1) - xDist;
        currentPk(2) = currentTileY(1) - yDist;
        storeAllPeaks{1,2} = currentPk;
        
        [~,~,peakStore] = tileUp(currentTileX,currentTileY,changeStore,numReps,xDist,yDist);
        storeAllPeaks{2,2} = peakStore;
        
        [~,~,peakStore] = tileDown(currentTileX,currentTileY,changeStore,numReps,xDist,yDist);
        storeAllPeaks{3,2} = peakStore;
        
        %% RIGHT
        currCol = size(storeAllPeaks,2);
        for iCol = 1:numReps
            newTileX = currentTileX + xChangeR;
            newTileY = currentTileY + yChangeR;
            
            plot(newTileX,newTileY,'b-')
            currentPk(1) = newTileX(1) - xDist;
            currentPk(2) = newTileY(1) - yDist;
            plot(currentPk(1),currentPk(2),'bo')
            storeAllPeaks{1,currCol+iCol} = currentPk;
            
            [~,~,peakStore] = tileUp(newTileX,newTileY,changeStore,numReps,xDist,yDist);
            storeAllPeaks{2,currCol+iCol} = peakStore;
            
            [~,~,peakStore] = tileDown(newTileX,newTileY,changeStore,numReps,xDist,yDist);
            storeAllPeaks{3,currCol+iCol} = peakStore;
            
            currentTileX = newTileX;
            currentTileY = newTileY;
        end
        
        %% LEFT
        currCol = size(storeAllPeaks,2);
        
        currentTileX = leftTile(:,1);
        currentTileY = leftTile(:,2);
        currentPk(1) = currentTileX(1) - xDist;
        currentPk(2) = currentTileY(1) - yDist;
        storeAllPeaks{1,currCol+1} = currentPk;
        
        [~,~,peakStore] = tileUp(currentTileX,currentTileY,changeStore,numReps,xDist,yDist);
        storeAllPeaks{2,currCol+1} = peakStore;
        
        [~,~,peakStore] = tileDown(currentTileX,currentTileY,changeStore,numReps,xDist,yDist);
        storeAllPeaks{3,currCol+1} = peakStore;
        
        %%
        currCol = size(storeAllPeaks,2);
        
        for iCol = 1:numReps
            newTileX = currentTileX + xChangeL;
            newTileY = currentTileY + yChangeL;
            
            plot(newTileX,newTileY,'b-')
            currentPk(1) = newTileX(1) - xDist;
            currentPk(2) = newTileY(1) - yDist;
            plot(currentPk(1),currentPk(2),'bo')
            storeAllPeaks{1,currCol+iCol} = currentPk;
            
            [~,~,peakStore] = tileUp(newTileX,newTileY,changeStore,numReps,xDist,yDist);
            storeAllPeaks{2,currCol+iCol} = peakStore;
            
            [~,~,peakStore] = tileDown(newTileX,newTileY,changeStore,numReps,xDist,yDist);
            storeAllPeaks{3,currCol+iCol} = peakStore;
            
            currentTileX = newTileX;
            currentTileY = newTileY;
        end
        
        %%
        hold on
        xline(40)
        yline(40)
        
        xline(0)
        yline(0)
        
        xline(80)
        yline(80)
        
        axis equal
        %     close(gcf)
        
        %% collapse pks into a column and plot them
        compilePks = [];
        figure('position',[1,41,1280,683]);
        for iRow = 1:size(storeAllPeaks,1)
            for iCol = 1:size(storeAllPeaks,2)
                if ~isempty(storeAllPeaks{iRow,iCol})
                    temp = storeAllPeaks{iRow,iCol};
                    plot(temp(:,1),temp(:,2),'o','color',colorMap(iMod,:))
                    hold on
                    compilePks = [compilePks; temp];
                end
            end
        end
        compilePks(end,:) = pkStore(iMod,:);
        hold on
        plot(pkStore(iMod,1),pkStore(iMod,2),'o','color',colorMap(iMod,:))
        boxCoords = [0 0; 0 80; 80 80; 80 0; 0 0];
        plot(boxCoords(:,1),boxCoords(:,2),'k-')
        axis equal
        xline(40)
        yline(40)
        
        modPks{1,iMod} = compilePks;
    end
    
    close all
    
    %% shift
    shiftStore = cell(1);
    finalHexPts = cell(1);
    finalModPts = cell(1);
    finalPtPairs = cell(1);
    distStore = nan(numRuns,3);
    
    % shift each module
    for iShift = 1:numRuns
        
        for iMod = 1:numMods
            temp = modPks{1,iMod};
            
            shiftedX = temp(:,1) + randi([-80 80],1);
            shiftedY = temp(:,2) + randi([-80 80],1);
            %       plot(shiftedX,shiftedY,'go')
            
            shiftStore{iShift,iMod} = [shiftedX shiftedY];
        end
        
        %% find the ones in the polygon for each module
        figure('position',[1,41,1280,683]);
        allPkStore = cell(1);
        flagStore = nan(1,numMods);
        
        for iMod = 1:numMods
            temp = shiftStore{iShift,iMod};
            modIdx = inpolygon(temp(:,1),temp(:,2),ccCoordStore(:,1,iMod),ccCoordStore(:,2,iMod));
            
            ptsToKeep = temp(modIdx,:);
            plot(ptsToKeep(:,1),ptsToKeep(:,2),'o','color',colorMap(iMod,:))
            hold on
            
            plot(ccCoordStore(:,1,iMod),ccCoordStore(:,2,iMod),'-','color',colorMap(iMod,:))
            allPkStore{1,iMod} = ptsToKeep;
            if isempty(ptsToKeep)
                flagStore(1,iMod) = 0;
            end
        end
        plot(boxCoords(:,1),boxCoords(:,2),'k-')
        xlim([-5 85])
        ylim([-5 85])
        axis square
        
%         close(gcf)
        
        finalHexPts{iShift,1} = allPkStore{1,1};
        finalHexPts{iShift,2} = allPkStore{1,2};
        finalHexPts{iShift,3} = allPkStore{1,3};
        
        %%
        if ~sum(flagStore == 0)
            
            if numMods == 2
                modComps = [1 2];
            else
                modComps = [1 2; 1 3; 2 3];
            end
            
            pltCnt = 1;
            for iComp = 1:size(modComps,1)
                for iRun = 1:2
                    % P1, P2: 1x2 vectors, coordinates inside hexagon
                    % V: Nx2 array, vertices of regular hexagon (used to compute tiling)
                    refMod = modComps(iComp,iRun);
                    
                    V = ccCoordStore(:,:,refMod);
                    P2 = allPkStore{1,refMod}; 
                    
                    % Step 1: Compute center of hex
                    center = [40 40];
                    
                    % Step 2: Get two edge directions from adjacent edges
                    dir1 = mean([V(1,:); V(2,:)]) - center;
                    dir2 = mean([V(6,:); V(1,:)]) - center;
                    
                    % Step 3: Construct lattice vectors
                    b1 = 2 * dir1;
                    b2 = 2 * dir2;
                    
                    % Step 4: Try wrapped versions of P2
                    wrapOffsets = [-1 -1; -1 0; -1 1;
                        0 -1;  0 0;  0 1;
                        1 -1;  1 0;  1 1];
                    
                    dists = zeros(size(wrapOffsets,1), 1);
                    ptStore = nan(size(wrapOffsets,1),2);
                    
                    for k = 1:size(wrapOffsets,1)
                        offset = wrapOffsets(k,1) * b1 + wrapOffsets(k,2) * b2;
                        P2_wrapped = P2 + offset;
                        ptStore(k,:) = P2_wrapped;
                    end
                    finalPtPairs{iComp,iRun} = ptStore;
                    
                    if PLOT == true
                        if iComp == 1 && iRun == 1
                            figure;
                            pltCnt = 1;
                        end
                        
                        subplot(1,size(modComps,1),pltCnt)
                        plot(V(:,1),V(:,2),'-','color',colorMap(refMod,:))
                        hold on
                        plot(ccCoordStore(:,1,refMod),ccCoordStore(:,2,refMod),'-','color',colorMap(refMod,:))
                        xline(40)
                        yline(40)
                    
                        plot(P2(1,1),P2(1,2),'o','color',colorMap(refMod,:))
                        plot(ptStore(:,1),ptStore(:,2),'o','color',colorMap(refMod,:))
                        
                        xlim([10 70])
                        ylim([10 70])
                        axis square
                    end
                end
                pltCnt = pltCnt + 1;
            end
            
            %%
            pltCnt = pltCnt - size(modComps,1);
            
            
            for iComp = 1:size(modComps,1)
                
                V = ccCoordStore(:,:,modComps(iComp,1));
                P = allPkStore{1,modComps(iComp,1)};

                [minDist, closestEdgeIdx, closestPoint, projPoints, dists] = closestEdgeToPoint(V,P);
                distToEdge = minDist*boxSize/szmap;
                
                distThresh = spacingAllExpt(iExp,modComps(iComp,1))*0.1;
                
                if distToEdge < distThresh
                    tempPts1 = finalPtPairs{iComp,1};
                    
                    distToPoly = nan(size(tempPts1,1),1);
                    for iPt = 1:size(tempPts1,1)
                        P = tempPts1(iPt,:);
                        [minDist, closestEdgeIdx, closestPoint, projPoints, dists] = closestEdgeToPoint(V,P);
                        distToPoly(iPt,1) = minDist*boxSize/szmap;
                    end
                    
                    inds = distToPoly < distThresh;
                    modPts1 = tempPts1(inds,:);
                else
                    modPts1 = P;
                end
                
                modPts2 = finalPtPairs{iComp,2};
                
                ptDist = nan(size(modPts1,1),size(modPts2,1));
                for iPt = 1:size(modPts1,1)
                    for jPt = 1:size(modPts2,1)
                        ptDist(iPt,jPt) = pdist([modPts1(iPt,:); modPts2(jPt,:)])*boxSize/szmap;
                    end
                end
                
                minVal = nanmin(nanmin(ptDist));
                [r,c] = find(ptDist == minVal);
                finalModPts{iShift,iComp} = [modPts1(r,:); modPts2(c,:)];
                
                if PLOT == true
                    subplot(1,size(modComps,1),pltCnt)
                    hold on
                    plot(modPts1(:,1),modPts1(:,2),'*')
                    plot(modPts2(:,1),modPts2(:,2),'*')
                    
                    plot([modPts1(r,1) modPts2(c,1)],[modPts1(r,2) modPts2(c,2)],'k-')
                    pltCnt = pltCnt + 1;
                end
                
                distStore(iShift,iComp) = minVal;
            end
        end
    end
    finalDistStore{iExp,1} = distStore;
    finalModPtStore{iExp,1} = finalModPts;
end

