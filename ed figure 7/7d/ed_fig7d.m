load expList_MEC_26718_rpt
load rotCCstore_26718_day1
numMods = size(rotCCstore,4);
mxStoreDay1 = mxStore;

%%
load colorMapLight
load colorMap 

% figure('position',[93,235,202,717]); one column
figure('position',[99,263,560,631]);
hold on
szmap = 40;
boxSize = 150;

cnt = 1;

% numMods = 3;

for iMod = 1:numMods
    subplot(3,3,cnt)
    
    if size(rotCCstore,3) == 3
        cc = rotCCstore(:,:,iMod);
    else
        cc = rotCCstore(:,:,mindStore(iMod),iMod);
    end
    
    colorMapBRK(cc);
    hold on
    plot([szmap szmap],[0 szmap*2],'k-','linew',1)
    plot([0 szmap*2],[szmap szmap],'k-','linew',1)
    
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
    
    pkStore(iMod,1,1) = x;
    pkStore(iMod,2,1) = y;
    
    plot(pkStore(iMod,1,1),pkStore(iMod,2,1),'k.','markersize',15)
    plot(ccCoordStore(:,1,iMod),ccCoordStore(:,2,iMod),'k-','linewidth',2)
    
    cnt = cnt + 3;
end

%%
colorblindMap
col = [3 6 9];

for iMod = 1:numMods
    
    subplot(3,3,col(iMod))

    plot(pkStore(iMod,1,1),pkStore(iMod,2,1),'.','markersize',25,'color',colorMapLight(iMod,:))
    hold on
    plot(ccCoordStore(:,1,iMod),ccCoordStore(:,2,iMod),'-','linewidth',2,'color',colorMapLight(iMod,:))
    set(gca,'ydir','rev')
    axis square
    xlim([15 65])
    ylim([15 65])
    set(gca,'xtick','','ytick','')
end

%%
cnt = 2;
load rotCCstore_26718_day2.mat
mxStoreDay2 = mxStore;

%%
for iMod = 1:numMods
    subplot(3,3,cnt)
    
    if size(rotCCstore,3) == 3
        cc = rotCCstore(:,:,iMod);
    else
        cc = rotCCstore(:,:,mindStore(iMod),iMod);
    end
    
    colorMapBRK(cc);
    hold on
    plot([szmap szmap],[0 szmap*2],'k-','linew',1)
    plot([0 szmap*2],[szmap szmap],'k-','linew',1)
    
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
    
    pkStore(iMod,1,2) = x;
    pkStore(iMod,2,2) = y;
    
    plot(pkStore(iMod,1,2),pkStore(iMod,2,2),'k.','markersize',15)
    plot(ccCoordStore(:,1,iMod),ccCoordStore(:,2,iMod),'k-','linewidth',2)
    
    cnt = cnt + 3;
end

%%
col = [3 6 9];

for iMod = 1:numMods
    
    subplot(3,3,col(iMod))

    plot(pkStore(iMod,1,2),pkStore(iMod,2,2),'.','markersize',25,'color',colorMap(iMod,:))
    hold on
    plot(ccCoordStore(:,1,iMod),ccCoordStore(:,2,iMod),'-','linewidth',2,'color',colorMap(iMod,:))
    set(gca,'ydir','rev')
    axis square
    xlim([15 65])
    ylim([15 65])
    set(gca,'xtick','','ytick','')
    
    dist(iMod,1) = pdist([pkStore(iMod,:,1); pkStore(iMod,:,2)])*boxSize/szmap;
    title(sprintf('Distance = %.2f cm',dist(iMod,1)))
    
end

%%
figure('position',[680,678,262,420]);
for iMod = 1:numMods
    subplot(numMods,1,iMod)
    plot(1:120,mxStoreDay1(iMod,:),'-','color',colorMapLight(iMod,:))
    hold on
    plot(1:120,mxStoreDay2(iMod,:),'-','color',colorMap(iMod,:))
ylabel('Correlation')

end    
box on
xlabel('Rotation (deg)')
set(gca,'xtick',0:20:120,'xticklabels',0:60:360)
% ylim([-0.005 0.305])
set(gca,'ytick',0:0.1:1)
