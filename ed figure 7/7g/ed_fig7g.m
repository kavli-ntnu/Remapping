load rotCCstore_26821_probe1
mxStoreP1 = mxStore;
maxCorrP1 = maxCorr;

%% these are added so that both probes are compared to combined mindStore 
mindStore(1) = 3;
mindStore(2) = 5;
mindStore(3) = 5;

%%
load colorMapLight
load colorMap

% figure('position',[93,235,202,717]); one column
figure('position',[99,263,560,631]);
hold on
szmap = 40;
boxSize = 150;

cnt = 1;

numMods = 3;

for iMod = 1:numMods
    subplot(3,1,cnt)
    
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
    
    cnt = cnt + 1;
end

ccCoordStoreP1 = ccCoordStore;
pkStoreP1 = pkStore;

%%
load rotCCstore_26821_probe2
mxStoreP2 = mxStore;
maxCorrP2 = maxCorr;

%%
mindStore(1) = 3;
mindStore(2) = 5;
mindStore(3) = 5;

%%
% figure('position',[93,235,202,717]); one column
figure('position',[99,263,560,631]);
hold on
szmap = 40;
boxSize = 150;

cnt = 1;

numMods = 3;

for iMod = 1:numMods
    subplot(3,1,cnt)
    
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
    
    cnt = cnt + 1;
end

ccCoordStoreP2 = ccCoordStore;
pkStoreP2 = pkStore;

%%
load rotCCstore_26821_bothProbes

%%
% figure('position',[93,235,202,717]); one column
figure('position',[99,263,560,631]);
hold on
szmap = 40;
boxSize = 150;

cnt = 1;

numMods = 3;

for iMod = 1:numMods
    subplot(3,1,cnt)
    
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
    
    cnt = cnt + 1;
end

%%
figure('position',[99,263,560,631]);

for iMod = 1:numMods
    
    subplot(3,1,iMod)
 
    plot(pkStore(iMod,1,1),pkStore(iMod,2,1),'.','markersize',25,'color',[0.5 0.5 0.5])
    hold on
    plot(ccCoordStore(:,1,iMod),ccCoordStore(:,2,iMod),'-','linewidth',2,'color',[0.5 0.5 0.5])
    
    
    plot(pkStoreP1(iMod,1,1),pkStoreP1(iMod,2,1),'.','markersize',25,'color',colorMapLight(iMod,:))
    hold on
    plot(ccCoordStoreP1(:,1,iMod),ccCoordStoreP1(:,2,iMod),'-','linewidth',2,'color',colorMapLight(iMod,:))
    
    plot(pkStoreP2(iMod,1,1),pkStoreP2(iMod,2,1),'.','markersize',25,'color',colorMap(iMod,:))
    hold on
    plot(ccCoordStoreP2(:,1,iMod),ccCoordStoreP2(:,2,iMod),'-','linewidth',2,'color',colorMap(iMod,:))
    
    set(gca,'ydir','rev')
    axis square
    xlim([15 65])
    ylim([15 65])
    set(gca,'xtick','','ytick','')
    yline(40)
    xline(40)
end

%%
plt = [1 2 3];

figure('position',[496,558,1435,340]);
for iMod = 1:numMods
    subplot(1,3,plt(iMod))
    plot(1:120,mxStoreP1(iMod,:),'-','color',colorMapLight(iMod,:),'linewidth',2)
    hold on
    set(gca,'xtick',0:20:120,'xticklabels',0:60:360)
    set(gca,'ytick',0:0.1:1)
    ylabel('Correlation')
    yMax(1,iMod) = maxCorr(iMod);
end


for iMod = 1:numMods
    subplot(1,3,plt(iMod))
    plot(1:120,mxStoreP2(iMod,:),'-','color',colorMap(iMod,:),'linewidth',2)
    hold on
    set(gca,'xtick',0:20:120,'xticklabels',0:60:360)
    set(gca,'ytick',0:0.1:1)
    yMax(2,iMod) = maxCorr(iMod);
end

for iMod = 1:numMods
    subplot(1,3,plt(iMod))
    plot(1:120,mxStore(iMod,:),'-','color',[0.5 0.5 0.5],'linewidth',2)
    hold on
    set(gca,'xtick',0:20:120,'xticklabels',0:60:360)
    set(gca,'ytick',0:0.1:1)
    ylabel('Correlation')
    yMax(3,iMod) = maxCorr(iMod);
end

subplot(131)
ylim([0 nanmax(yMax(:,1))+0.05])
subplot(132)
ylim([0 nanmax(yMax(:,2))+0.05])
subplot(133)
ylim([0 nanmax(yMax(:,3))+0.05])
