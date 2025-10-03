load expList_26718
load rotCCstore_26718
numMods = size(rotCCstore,4);
load rotCCstore_singleAB_26718

boxSize = expList_MEC{1,7}(1);
szmap = 40;

%%
pkStore = nan(numMods,2);
maskStore = nan(size(rotCCstore,1),size(rotCCstore,1),numMods);
ccCoordStore = nan(7,2,numMods);
spacingStore = nan(1,numMods);
orientStore = nan(numMods,3);

for iMod = 1:numMods
    cc = rotCCstore(:,:,mindStore(iMod),iMod);
    
    cc(isnan(cc)) = 0;
    ac = xcorr2(cc);
    [gScore,gStats,~,COM] = analyses.gridnessScore(ac);
    origin = round(size(ac,1) / 2);
    
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
    %         colorMapBRK(rotCCstore(:,:,iMod));
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
    
    close(gcf);
end

%% plot points for individual modules
for iMod = 1:numMods
    
    ccStoreAll = rotCCstore_singleAB_26718{1,iMod};

    figure;
    hold on
    
    fieldStore = [];
    for iCell = 1:size(ccStoreAll,4)
        
        if ~isempty(ccStoreAll(:,:,mindStore(iMod),iCell))
            cc = ccStoreAll(:,:,mindStore(iMod),iCell);
            [f1,f2] = analyses.placefield(cc,'threshold',0.1,'minPeak',0.1);
            
            if ~isempty(f2)
                fieldTemp = nan(size(f2,2),2);
                for iField = 1:size(f2,2)
                    plot(f2(iField).peakX,f2(iField).peakY,'ko')
                    fieldTemp(iField,1) = f2(iField).peakX;
                    fieldTemp(iField,2) = f2(iField).peakY;
                end
                
                fieldStore = [fieldStore; fieldTemp];
            end
        end
    end
    
    szmap = 40;
    plot([szmap szmap],[0 szmap*2],'k-','linew',1)
    plot([0 szmap*2],[szmap szmap],'k-','linew',1)
    
    close(gcf)
    
    %%
    uniqueLoc = unique(fieldStore,'rows');
    
    numTimes = nan(size(uniqueLoc,1),1);
    for iLoc = 1:size(uniqueLoc,1)
        tempLoc = uniqueLoc(iLoc,:);
        
        numX = sum(fieldStore(:,1) == tempLoc(1) & fieldStore(:,2) == tempLoc(2));
        
        numTimes(iLoc,1) = numX;
    end
    
    %%
    figure;
    hold on
    
    uniqueColors = unique(numTimes);
    numSteps = size(uniqueColors,1);
    cmap = makeColorMap([1 1 1],[0 0 0],size(unique(numTimes),1));
    
    colorRow = nan(size(numTimes,1),1);
    for iRow = 1:size(numTimes,1)
        temp = numTimes(iRow,1);
        colorRow(iRow,1) = find(uniqueColors == temp);
    end
    
    for iLoc = 1:size(uniqueLoc,1)
        %         colorX = 1 - 0.09 * numTimes(iLoc);
        %         plot(uniqueLoc(iLoc,1),uniqueLoc(iLoc,2),'ko','markerFaceColor',[colorX colorX colorX])
        plot(uniqueLoc(iLoc,1),uniqueLoc(iLoc,2),'ko','markerFaceColor',cmap(colorRow(iLoc),:))
    end
    
    szmap = 40;
    plot([szmap szmap],[0 szmap*2],'k-','linew',1)
    plot([0 szmap*2],[szmap szmap],'k-','linew',1)
    axis square
    set(gca, 'ydir', 'rev')
    
    hold on
    plot(ccCoordStore(:,1,iMod),ccCoordStore(:,2,iMod),'k-','linewidth',2)
    xlim([0 80])
    ylim([0 80])
    
end
