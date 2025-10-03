
load expList_26821
load gridIDstore_26821
numMods = size(gridIDstore,2);
load rotCCstore_AB_26821

distStore = nan(size(expList_MEC,1),3);
pltCnt = 1;

for iExp = 1
    
    exptSelected = iExp;
    
    name = expList_MEC{exptSelected,3};
    boxSize = expList_MEC{exptSelected,7}(1);
    szmap = 40;
    
    %%
    pkStore = nan(numMods,2);
    maskStore = nan(size(rotCCstore,1),size(rotCCstore,1),numMods);
    ccCoordStore = nan(7,2,numMods);
    spacingStore = nan(1,numMods);
    orientStore = nan(numMods,3);
    
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
            prct = 0.1;
            
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
    
    newPg = 1:30:1000;
    if sum(pltCnt == newPg)
         figure('position',[-1919,-49,1920,955]);
         pltCnt = 1;
    end
    
    for iMod = 1:numMods
        subplot(5,6,pltCnt)
        %     colorMapBRK(rotCCstore(:,:,iMod));
        colorMapBRK(rotCCstore(:,:,mindStore(iMod),iMod));
        hold on
        %     plot(winners(iMod,1),winners(iMod,2),'k.','markersize',15)
        plot(pkStore(iMod,1),pkStore(iMod,2),'k.','markersize',15)
        xline(40)
        yline(40)
        %     xlim([15 65])
        %     ylim([15 65])
        plot(ccCoordStore(:,1,iMod),ccCoordStore(:,2,iMod),'-','linewidth',2,'color',colorMap(iMod,:))
        pltCnt = pltCnt + 1;
    end
    
    %%
    if numMods == 2
        modComps = [1 2];
    elseif numMods == 3
        modComps = [1 2; 1 3; 2 3];
    else
        modComps = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    end
    
    for iComp = 1:size(modComps,1)
        for iRun = 1:2
            % P1, P2: 1x2 vectors, coordinates inside hexagon
            % V: Nx2 array, vertices of regular hexagon (used to compute tiling)
            refMod = modComps(iComp,iRun);
            
            V = ccCoordStore(:,:,refMod);
            P2 = pkStore(refMod,:);
            
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
            
            subplot(5,6,pltCnt)
            plot(V(:,1),V(:,2),'-','color',colorMap(refMod,:))
            hold on
            plot(ccCoordStore(:,1,refMod),ccCoordStore(:,2,refMod),'-','color',colorMap(refMod,:))
            xline(40)
            yline(40)
            
            %         plot(P1(1,1),P1(1,2),'o','color',colorMap(testMod,:))
            plot(P2(1,1),P2(1,2),'o','color',colorMap(refMod,:))
            set(gca,'ydir','rev')
            
            plot(ptStore(:,1),ptStore(:,2),'o','color',colorMap(refMod,:))
            
            xlim([10 70])
            ylim([10 70])
            axis square
        end 
        pltCnt = pltCnt + 1;
    end
    
    %%
    pltCnt = pltCnt - size(modComps,1);
    
    finalModPts = cell(1);
    for iComp = 1:size(modComps,1)   
        
        V = ccCoordStore(:,:,modComps(iComp,1));
        P = pkStore(modComps(iComp,1),:);
        [minDist, closestEdgeIdx, closestPoint, projPoints, dists] = closestEdgeToPoint(V,P);
        distToEdge = minDist*boxSize/szmap;
        
        distThresh = spacingStore(iExp,modComps(iComp,1))*0.1;
        
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
        finalModPts{1,iComp} = [modPts1(r,:); modPts2(c,:)];
        
        subplot(5,6,pltCnt)
        hold on
        plot(modPts1(:,1),modPts1(:,2),'*')
        plot(modPts2(:,1),modPts2(:,2),'*')
        
        plot([modPts1(r,1) modPts2(c,1)],[modPts1(r,2) modPts2(c,2)],'k-')
        pltCnt = pltCnt + 1;
%         axis off

        distStore(iExp,iComp) = minVal;
    end
end

distStore