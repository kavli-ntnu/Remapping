load expList_MEC_28367_FN
load rotCCstore_28367_0923
mxStore1 = mxStore;

iExp = 1; exptSelected = iExp;

numMods = size(rotCCstore,4);
selectedMods = 1:numMods;
boxSize = expList_MEC{1,7}(1);

load colorMap
load colorMapLight

%% loop through each module to get pts of ac
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
    
    %     szmap = 39;
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
    
    close(gcf)
end

% close all

allPkStore23 = allPkStore;
modSpacing23 = modSpacing;

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
%     origin = round(size(ac,1) / 2);
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
    close(gcf)
end

pkStore23 = pkStore;
ccStore23 = ccCoordStore;

%%
figure;
pltCnt = 1;
for iMod = 1:3
    subplot(3,3,pltCnt)
    colorMapBRK(rotCCstore(:,:,mindStore(iMod),iMod));
    hold on
    plot(ccStore23(:,1,iMod),ccStore23(:,2,iMod),'-','color',colorMap(iMod,:),'linewidth',2)
    plot(pkStore23(iMod,1),pkStore23(iMod,2),'kx','markersize',10)
    xline(40)
    yline(40)
    title(sprintf('%.2f',maxCorr(1,iMod)))
    pltCnt = pltCnt+3;
end

%%
iExp = 2;
load rotCCstore_28367_0926
% mindStore = finalMindStore;
mxStore2 = mxStore;

%% loop through each module to get pts of ac
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
    
    %     szmap = 39;
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
    
    
    insideStore{1,iMod} = [vxInside vyInside];
    polyStore(:,1:2,iMod) = [acX_converted acY_converted];
    allPkStore{1,iMod} = [vx vy];
    
            close(gcf)
end

% close all

allPkStore26 = allPkStore;
modSpacing26 = modSpacing;

%
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
%     origin = round(size(ac,1) / 2);
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
    close(gcf)
end

pkStore26 = pkStore;
ccStore26 = ccCoordStore;

%%
pltCnt = 2;
for iMod = 1:3
    subplot(3,3,pltCnt)
    colorMapBRK(rotCCstore(:,:,mindStore(iMod),iMod));
    hold on
    plot(ccStore26(:,1,iMod),ccStore26(:,2,iMod),'-','color',colorMap(iMod,:),'linewidth',2)
    plot(pkStore26(iMod,1),pkStore26(iMod,2),'kx','markersize',10)
    xline(40)
    yline(40)
    title(sprintf('%.2f',maxCorr(1,iMod)))
    pltCnt = pltCnt+3;
end

%%
iExp = 3;
load rotCCstore_28367_0928
% mindStore = finalMindStore;
mxStore3 = mxStore;

%% loop through each module to get pts of ac
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
    
    %     szmap = 39;
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
    close(gcf)
end

% close all

allPkStore28 = allPkStore;
modSpacing28 = modSpacing;

%
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
%     origin = round(size(ac,1) / 2);
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
    close(gcf)
end

pkStore28 = pkStore;
ccStore28 = ccCoordStore;

%%
pltCnt = 3;
for iMod = 1:3
    subplot(3,3,pltCnt)
    colorMapBRK(rotCCstore(:,:,mindStore(iMod),iMod));
    hold on
    plot(ccStore28(:,1,iMod),ccStore28(:,2,iMod),'-','color',colorMap(iMod,:),'linewidth',2)
    plot(pkStore28(iMod,1),pkStore28(iMod,2),'kx','markersize',10)
    xline(40)
    yline(40)
    title(sprintf('%.2f',maxCorr(1,iMod)))
    pltCnt = pltCnt+3;
end


%%
figure;
for iMod = selectedMods
    subplot(numMods,1,iMod)
    plot(ccStore28(:,1,iMod),ccStore28(:,2,iMod),'-','linewidth',2,'color',colorMap(iMod,:))
    hold on
    plot(ccStore26(:,1,iMod),ccStore26(:,2,iMod),'-','linewidth',2,'color',colorMapLight(iMod,:))
    plot(ccStore23(:,1,iMod),ccStore23(:,2,iMod),'-','linewidth',2,'color',[0.5 0.5 0.5])   
    
    pks28 = allPkStore28{1,iMod};
    plot(pks28(:,1),pks28(:,2),'o','color',colorMap(iMod,:),'linewidth',1.5,'markersize',5)
    
    pks26 = allPkStore26{1,iMod};
    plot(pks26(:,1),pks26(:,2),'*','color',colorMapLight(iMod,:),'linewidth',1.5,'markersize',7)
    
    pks23 = allPkStore23{1,iMod};
    plot(pks23(:,1),pks23(:,2),'.','color',[0.5 0.5 0.5],'linewidth',1,'markersize',10)
    
    axis square
    xline(40)
    yline(40)
    set(gca,'ydir','rev')
    ylim([0 80])
    xlim([0 80])
    set(gca,'xtick','','ytick','')
    
    if iMod == 1
        xlim([31 49])
        ylim([31 49]) % m1 zoom
    elseif iMod == 2
        xlim([25 55]) % m2 zoom
        ylim([25 55])
    elseif iMod == 3
        xlim([20 60]) % m3 zoom
        ylim([20 60])
    end

end

%%
dM1(1) = pdist([pkStore23(1,:); pkStore26(1,:)])*3.75;
dM1(2) = pdist([pkStore23(1,:); pkStore28(1,:)])*3.75;
dM1(3) = pdist([pkStore26(1,:); pkStore28(1,:)])*3.75

nanmean(dM1)

dM2(1) = pdist([pkStore23(2,:); pkStore26(2,:)])*3.75;
dM2(2) = pdist([pkStore23(2,:); pkStore28(2,:)])*3.75;
dM2(3) = pdist([pkStore26(2,:); pkStore28(2,:)])*3.75

nanmean(dM2)

dM3(1) = pdist([pkStore23(3,:); 52 47])*3.75;
dM3(2) = pdist([pkStore23(3,:); pkStore28(3,:)])*3.75;
dM3(3) = pdist([52 47; pkStore28(3,:)])*3.75

nanmean(dM3)

%%
mxStore = mxStore1;

% smooth mod 1
for iCol = 1:120
    if iCol == 42 
        mxStore(1,iCol) = (mxStore(1,iCol-2)+mxStore(1,iCol-1)) / 2;
    elseif iCol == 43
        mxStore(1,iCol) = (mxStore(1,iCol+1)+mxStore(1,iCol+2)) / 2;
    elseif mxStore(1,iCol) < 0.1
        mxStore(1,iCol) = (mxStore(1,iCol-1)+mxStore(1,iCol+1)) / 2;
    end
end

figure('position',[-688,22,281,741]);
for iMod = 1:numMods
    subplot(3,1,iMod)
    plot(1:120,mxStore(iMod,:),'color',[0.5 0.5 0.5],'linewidth',2)
    set(gca,'xtick',0:20:120,'xticklabels',0:60:360)
    set(gca,'ytick',0:0.1:0.4)
    ylim([0 0.4])
end

%
mxStore = mxStore2;

for iMod = 1:numMods
    subplot(3,1,iMod)
    hold on
    plot(1:120,mxStore(iMod,:),'color',colorMapLight(iMod,:),'linewidth',2)
    set(gca,'xtick',0:20:120,'xticklabels',0:60:360)
    set(gca,'ytick',0:0.1:round(maxCorr(iMod),1))
    set(gca,'ytick',0:0.1:0.4)
end

%
mxStore = mxStore3;

% smooth mod 1
for iCol = 1:120
    if iCol == 42 
        mxStore(1,iCol) = (mxStore(1,iCol-2)+mxStore(1,iCol-1)) / 2;
    elseif iCol == 43
        mxStore(1,iCol) = (mxStore(1,iCol+1)+mxStore(1,iCol+2)) / 2;
    elseif mxStore(1,iCol) < 0.1
        mxStore(1,iCol) = (mxStore(1,iCol-1)+mxStore(1,iCol+1)) / 2;
    end
end

for iMod = 1:numMods
    subplot(3,1,iMod)
    hold on
    plot(1:120,mxStore(iMod,:),'color',colorMap(iMod,:),'linewidth',2)
    set(gca,'xtick',0:20:120,'xticklabels',0:60:360)
    set(gca,'ytick',0:0.1:0.4)
    ylim([0 0.4])
    ylabel('PV correlation')
end

