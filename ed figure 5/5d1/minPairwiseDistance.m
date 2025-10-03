load expList_26718
load rotCCstore_26718
numMods = size(rotCCstore,4);

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
            orientStore(iMod,1:3) = (orientation(1:3,1))';
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

    %% set variables
    selectedMods = 1:numMods;
    numMods = size(selectedMods,2);
    
    % loop through each module to get pts of ac
    szmap = 40;
    modSpacing = nan(size(selectedMods,2),1);
    insideStore = cell(1);
    allPkStore = cell(1);
    shiftedPks = cell(1);
    % polyStore = nan(6,2,numMods);
    
    for iMod = selectedMods
        cc = rotCCstore(:,:,mindStore(iMod),iMod);
        %         cc = rotCCstore(:,:,iMod); % expList_MEC_other coco
        cc(isnan(cc)) = 0;
        ac = xcorr2(cc);
        
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
        
        [~,stats] = analyses.gridnessScore(ac);
        modSpacing(iMod,1) = nanmean(stats.spacing) * boxSize/szmap;
        
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

        % if you want to check that ac points are used to create
        % polygon that is mapped onto cc
        figure;
        subplot(121)
        colorMapBRK(ac);
        hold on
        plot(acX,acY,'w.','markersize',10)
        axis on
        
        subplot(122)
        colorMapBRK(cc);
        hold on
        plot([szmap szmap],[0 szmap*2],'w-','linew',1)
        plot([0 szmap*2],[szmap szmap],'w-','linew',1)
        axis on
        
        acX_converted = acX - floor(max(size(cc))/2);
        acY_converted = acY - floor(max(size(cc))/2);
        plot(acX_converted,acY_converted,'w.','markersize',15)
        
        % find all peaks in ensemble CC
        cc(cc < 0.02) = eps;
        rm = imregionalmax(cc);
        rp = regionprops(rm);
        numPeaks = size(rp,1);
        prct = 0.2;
        
        % coordinates of peaks from cc
        vx = []; vy = [];
        for i = 1:length(rp)
            vx = [vx; rp(i).Centroid(1)];
            vy = [vy; rp(i).Centroid(2)];
        end
        
        % plot the peaks
        hold on
        plot(vx,vy,'wo')
        
        acX_converted = acX_converted(2:end) - floor(max(size(cc))/2);
        acY_converted = acY_converted(2:end) - floor(max(size(cc))/2);
        
        % sort points of ac sequentially
        [~, sind] = sort(atan2(acY_converted,acX_converted),'ascend');
        acX_converted = acX_converted(sind);
        acY_converted = acY_converted(sind);
        
        acX_converted = [acX_converted; acX_converted(1)];
        acY_converted = [acY_converted; acY_converted(1)];
        
        acX_converted = acX_converted + floor(max(size(cc))/2);
        acY_converted = acY_converted + floor(max(size(cc))/2);
        
        %         figure; plot(acX_converted,acY_converted)
        
        % find the peaks inside the polygon defined by acX and acY
        polyIdx = inpolygon(vx,vy,acX_converted,acY_converted);
        sum(polyIdx);
        % get the pts of vx and vy that are inside the polygon
        vxInside = vx(polyIdx);
        vyInside = vy(polyIdx);
        
        %     plot(vxInside,vyInside,'wo','markersize',15)
        %     saveas(gcf,fullfile(outFolder,'ac_cc.jpg'))
        %     close(gcf)
        
        insideStore{1,iMod} = [vxInside vyInside];
        polyStore(:,1:2,iMod) = [acX_converted acY_converted];
        allPkStore{1,iMod} = [vx vy];
        
        %     saveas(gcf,fullfile(exptFolder,sprintf('mod%d_ac_cc.jpg',iMod)))
        %     close(gcf)  
    end
    
    %%
    figure;
    for iMod = 1:numMods
        
        modPks = allPkStore{1,iMod};
        
        plot(modPks(:,1),modPks(:,2),'o','linewidth',1,'color',colorMap(iMod,:))
        hold on
        plot(ccCoordStore(:,1,iMod),ccCoordStore(:,2,iMod),'-','linewidth',2,'color',colorMap(iMod,:))
        
        axis square
        xlim([0 80])
        ylim([0 80])
        xline(40)
        yline(40)
        set(gca,'ydir','rev')
        set(gca,'xticklabels','','yticklabels','')
    end
       
    %%
    minPts = [];
    if numMods == 2
        averageDist = nan(size(allPkStore{1,selectedMods(end)}(:,1),1),size(allPkStore{1,selectedMods(end-1)}(:,1),1));
        distFromCenter = nan(size(allPkStore{1,selectedMods(end)}(:,1),1),size(allPkStore{1,selectedMods(end-1)}(:,1),1));
        
        for i = 1:size(allPkStore{1,selectedMods(end-1)}(:,1),1)
            for j = 1:size(allPkStore{1,selectedMods(end)}(:,1),1)
                averageDist(j,i) = pdist([allPkStore{1,selectedMods(end-1)}(i,1) allPkStore{1,selectedMods(end-1)}(i,2); allPkStore{1,selectedMods(end)}(j,1) allPkStore{1,selectedMods(end)}(j,2)]) * boxSize/szmap;
                
                
                dist1 = pdist([allPkStore{1,selectedMods(end-1)}(i,1) allPkStore{1,selectedMods(end-1)}(i,2); 39 39])* boxSize/szmap;
                dist2 = pdist([allPkStore{1,selectedMods(end)}(j,1) allPkStore{1,selectedMods(end)}(j,2); 39 39])* boxSize/szmap;
                distFromCenter(j,i) = nanmean([dist1 dist2]);
            end
        end
        
        minVal = nanmin(nanmin((averageDist)));
        fIdx = find(averageDist == minVal);
        [r,c] = ind2sub(size(averageDist),fIdx);
        
        minPts(1,:) = allPkStore{1,1}(c(1),:);
        minPts(2,:) = allPkStore{1,2}(r(1),:);
        
        for iMod = 1:numMods
            plot(minPts(iMod,1),minPts(iMod,2),'o','color',colorMap(iMod,:),'linewidth',2)
        end
        
    elseif numMods == 3
        averageDist = nan(size(allPkStore{1,selectedMods(end)}(:,1),1),size(allPkStore{1,selectedMods(end-1)}(:,1),1),size(allPkStore{1,selectedMods(1)}(:,1),1));
        distFromCenter = nan(size(allPkStore{1,selectedMods(end)}(:,1),1),size(allPkStore{1,selectedMods(end-1)}(:,1),1),size(allPkStore{1,selectedMods(1)}(:,1),1));
        
        for i = 1:size(allPkStore{1,selectedMods(1)}(:,1),1)
            for j = 1:size(allPkStore{1,selectedMods(end-1)}(:,1),1)
                for k = 1:size(allPkStore{1,selectedMods(end)}(:,1),1)
                    possDist1 = pdist([allPkStore{1,selectedMods(1)}(i,1) allPkStore{1,selectedMods(1)}(i,2); allPkStore{1,selectedMods(end-1)}(j,1) allPkStore{1,selectedMods(end-1)}(j,2)]) * boxSize/szmap;
                    possDist2 = pdist([allPkStore{1,selectedMods(1)}(i,1) allPkStore{1,selectedMods(1)}(i,2); allPkStore{1,selectedMods(end)}(k,1) allPkStore{1,selectedMods(end)}(k,2)]) * boxSize/szmap;
                    possDist3 = pdist([allPkStore{1,selectedMods(end-1)}(j,1) allPkStore{1,selectedMods(end-1)}(j,2); allPkStore{1,selectedMods(end)}(k,1) allPkStore{1,selectedMods(end)}(k,2)]) * boxSize/szmap;
                    averageDist(k,j,i) = nanmean([possDist1 possDist2 possDist3]);
                    
                    % edited correctly
                    dist1 = pdist([allPkStore{1,selectedMods(1)}(i,1) allPkStore{1,selectedMods(1)}(i,2); 39 39])* boxSize/szmap;
                    dist2 = pdist([allPkStore{1,selectedMods(end-1)}(j,1) allPkStore{1,selectedMods(end-1)}(j,2); 39 39])* boxSize/szmap;
                    dist3 = pdist([allPkStore{1,selectedMods(end)}(k,1) allPkStore{1,selectedMods(end)}(k,2); 39 39])* boxSize/szmap;
                    distFromCenter(k,j,i) = nanmean([dist1 dist2 dist3]);
                end
            end
        end
        
        minVal = nanmin(nanmin(nanmin(averageDist)));
        fIdx = find(averageDist == minVal);
        [r,c,d] = ind2sub(size(averageDist),fIdx);
        
        minPts(1,:) = allPkStore{1,1}(d,:);
        minPts(2,:) = allPkStore{1,2}(c,:);
        minPts(3,:) = allPkStore{1,3}(r,:);
        
        for iMod = 1:3
            plot(minPts(iMod,1),minPts(iMod,2),'o','color',colorMap(iMod,:),'linewidth',2)
        end
        
    elseif numMods == 4
        averageDist = nan(size(allPkStore{1,selectedMods(4)}(:,1),1),size(allPkStore{1,selectedMods(3)}(:,1),1),size(allPkStore{1,selectedMods(2)}(:,1),1),size(allPkStore{1,selectedMods(1)}(:,1),1));
        distFromCenter = nan(size(allPkStore{1,selectedMods(4)}(:,1),1),size(allPkStore{1,selectedMods(3)}(:,1),1),size(allPkStore{1,selectedMods(2)}(:,1),1),size(allPkStore{1,selectedMods(1)}(:,1),1));
        
        for i = 1:size(allPkStore{1,selectedMods(1)}(:,1),1)
            for j = 1:size(allPkStore{1,selectedMods(2)}(:,1),1)
                for k = 1:size(allPkStore{1,selectedMods(3)}(:,1),1)
                    for l = 1:size(allPkStore{1,selectedMods(4)}(:,1),1)
                        possDist1 = pdist([allPkStore{1,selectedMods(1)}(i,1) allPkStore{1,selectedMods(1)}(i,2); allPkStore{1,selectedMods(2)}(j,1) allPkStore{1,selectedMods(2)}(j,2)]) * boxSize/szmap;
                        possDist2 = pdist([allPkStore{1,selectedMods(1)}(i,1) allPkStore{1,selectedMods(1)}(i,2); allPkStore{1,selectedMods(3)}(k,1) allPkStore{1,selectedMods(3)}(k,2)]) * boxSize/szmap;
                        possDist3 = pdist([allPkStore{1,selectedMods(1)}(i,1) allPkStore{1,selectedMods(1)}(i,2); allPkStore{1,selectedMods(4)}(l,1) allPkStore{1,selectedMods(4)}(l,2)]) * boxSize/szmap;
                        
                        possDist4 = pdist([allPkStore{1,selectedMods(2)}(j,1) allPkStore{1,selectedMods(2)}(j,2); allPkStore{1,selectedMods(3)}(k,1) allPkStore{1,selectedMods(3)}(k,2)]) * boxSize/szmap;
                        possDist5 = pdist([allPkStore{1,selectedMods(2)}(j,1) allPkStore{1,selectedMods(2)}(j,2); allPkStore{1,selectedMods(4)}(l,1) allPkStore{1,selectedMods(4)}(l,2)]) * boxSize/szmap;
                        
                        possDist6 = pdist([allPkStore{1,selectedMods(3)}(k,1) allPkStore{1,selectedMods(3)}(k,2); allPkStore{1,selectedMods(4)}(l,1) allPkStore{1,selectedMods(4)}(l,2)]) * boxSize/szmap;
                        
                        averageDist(l,k,j,i) = nanmean([possDist1 possDist2 possDist3 possDist4 possDist5 possDist6]);
                        
                        dist1 = pdist([allPkStore{1,selectedMods(1)}(i,1) allPkStore{1,selectedMods(1)}(i,2); 39 39]) * boxSize/szmap;
                        dist2 = pdist([allPkStore{1,selectedMods(2)}(j,1) allPkStore{1,selectedMods(2)}(j,2); 39 39]) * boxSize/szmap;
                        dist3 = pdist([allPkStore{1,selectedMods(3)}(k,1) allPkStore{1,selectedMods(3)}(k,2); 39 39]) * boxSize/szmap;
                        dist4 = pdist([allPkStore{1,selectedMods(4)}(l,1) allPkStore{1,selectedMods(4)}(l,2); 39 39]) * boxSize/szmap;
                        
                        distFromCenter(l,k,j,i) = nanmean([dist1 dist2 dist3 dist4]);
                    end
                end
            end
        end
        
        minVal = nanmin(nanmin(nanmin(nanmin(averageDist))));
        fIdx = find(averageDist == minVal);
        [r,c,d,e] = ind2sub(size(averageDist),fIdx);
        
        minPts(1,:) = allPkStore{1,1}(e,:);
        minPts(2,:) = allPkStore{1,2}(d,:);
        minPts(3,:) = allPkStore{1,3}(c,:);
        minPts(4,:) = allPkStore{1,4}(r,:);
        
        for iMod = 1:numMods
            plot(minPts(iMod,1),minPts(iMod,2),'o','color',colorMap(iMod,:),'linewidth',2)
        end 
    end
    
    if sum(sum(minPts < 15))
        if sum(minPts(:,1) < 15)
            f15 = 15-nanmin(minPts(:,1))+3;
            xlim([15-f15 65+f15])
            ylim([15-f15 65+f15])
        elseif sum(minPts(:,2) < 15)
            f15 = 15-nanmin(minPts(:,2))+3;
            xlim([15-f15 65+f15])
            ylim([15-f15 65+f15])
        end
    else
        xlim([15 65])
        ylim([15 65])
    end
    title(sprintf('d = %.1f cm',minVal))
end

% minimum distance
minVal

minDistances = [];
if numMods == 2
    minDistances(1) = pdist([minPts(1,:); minPts(2,:)])*boxSize/szmap
elseif numMods == 3
    minDistances(1) = pdist([minPts(1,:); minPts(2,:)])*boxSize/szmap;
    minDistances(2) = pdist([minPts(1,:); minPts(3,:)])*boxSize/szmap;
    minDistances(3) = pdist([minPts(2,:); minPts(3,:)])*boxSize/szmap
else
    minDistances(1) = pdist([minPts(1,:); minPts(2,:)])*boxSize/szmap;
    minDistances(2) = pdist([minPts(1,:); minPts(3,:)])*boxSize/szmap;
    minDistances(3) = pdist([minPts(1,:); minPts(4,:)])*boxSize/szmap;
    minDistances(4) = pdist([minPts(2,:); minPts(3,:)])*boxSize/szmap;
    minDistances(5) = pdist([minPts(2,:); minPts(4,:)])*boxSize/szmap;
    minDistances(6) = pdist([minPts(3,:); minPts(4,:)])*boxSize/szmap
end
