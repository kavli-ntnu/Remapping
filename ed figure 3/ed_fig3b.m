
load expList_26821
load rotCCstore_AB_26821
numMods = size(rotCCstore,4);

boxSize = expList_MEC{1,7}(1);
szmap = round(size(rotCCstore,1)/2);

%%
pkStore = nan(numMods,2);
maskStore = nan(size(rotCCstore,1),size(rotCCstore,1),numMods);
ccCoordStore = nan(7,2,numMods);
spacingStore = nan(1,numMods);
orientStore = nan(numMods,3);

for iMod = 2 %1:numMods

    cc = rotCCstore(:,:,mindStore(iMod),iMod);
    
    cc(isnan(cc)) = 0;
    ac = xcorr2(cc);
    [gScore,gStats,~,COM] = analyses.gridnessScore(ac);
    origin = size(ac,1) / 2;

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
        if size(orientation,1) == 3
            hOrient = orientation + 30;
            orientStore(iMod,1:3) = (orientation)';
        elseif size(orientation,1) == 6
            o = orientation([1 3 5]);
            hOrient = o + 30;
            orientStore(iMod,1:3) = (o)';
        end
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
    colorMapBRK(ac,'ydir','normal');
    axis on
    
    hold on
    plot(COM(:,1),COM(:,2),'k.')
   
    subplot(1,3,2)
    colorMapBRK(rotCCstore(:,:,mindStore(iMod),iMod),'ydir','normal');
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
    colorMapBRK(xcmask,'ydir','normal');
    axis on
    hold on
    plot(x,y,'ko','markers',15)
    
    pkStore(iMod,1) = x;
    pkStore(iMod,2) = y;
    
    maskStore(:,:,iMod) = xcmask;
end

