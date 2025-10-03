
load expList_HP_29730
load sData_HP_29730
load placeCells_29730

iExp = 1;
sessionNums = expList_HP{iExp,5};
numSessions = size(sessionNums,2);
sessionComp = expList_HP{iExp,8};
numComps = size(sessionComp,1);
depthFilter = [];
smth = 1;

%% set box size and nDeg
boxSize = expList_HP{iExp,7}(1);
nBins = 40;
nDeg = 3;

%%
mxStore = nan(numComps,360/nDeg,2);
maxCorr = nan(numComps,2);
minCorr = nan(numComps,2);
rotAngle = nan(numComps,2);
mindStore = nan(numComps,2);
numStore = nan(numComps,2);
rotCCstore = nan(nBins*2-1,nBins*2-1,360/nDeg,numComps);
unitStore = cell(1);

for iComp = 1:size(sessionComp,1)
    
    % set units   
    units = placeCells';
    
    % make map stacks
    % rescale pos data to fit box and to match each other
    posT1 = sData{1,sessionComp(iComp,1)}.results.posT; % timestamps for X and Y
    posX1 = minions.rescaleData(sData{1,sessionComp(iComp,1)}.results.posX,0,boxSize); % positions X
    posY1 = minions.rescaleData(sData{1,sessionComp(iComp,1)}.results.posY,0,boxSize); % positions Y
    min1x = nanmin(posX1);
    max1x = nanmax(posX1);
    min1y = nanmin(posY1);
    max1y = nanmax(posY1);
    posT2 = sData{1,sessionComp(iComp,2)}.results.posT;
    posX2 = minions.rescaleData(sData{1,sessionComp(iComp,2)}.results.posX,min1x,max1x);
    posY2 = minions.rescaleData(sData{1,sessionComp(iComp,2)}.results.posY,min1y,max1y);
    
    cnt = 1;
    clear m1all m2all
    
    % maps are rescaled from 0 to 1 to ignore rate diffs and nans are filled in
    for iUnit = units
        
        if ~isempty(sData{1,sessionComp(iComp,1)}.results.spkPos{iUnit})
            spkT = sData{1,sessionComp(iComp,1)}.results.spkPos{iUnit}(:,1); % get spike positions for this cell in session 1
            map = analyses.map([posT1,posX1,posY1],spkT,'binWidth',boxSize/nBins,'smooth',smth); % plot position data with spike times
            m1 = minions.rescaleData(map.z,eps,1); % m1 = rescaled rate map for this cell in first session, avoid having 0s in map
            m1 = inpaint_nans(m1); % interpolate surrounding values to fill in nans
            
            if ~isempty(sData{1,sessionComp(iComp,2)}.results.spkPos{iUnit})
                spkT = sData{1,sessionComp(iComp,2)}.results.spkPos{iUnit}(:,1);
                map = analyses.map([posT2,posX2,posY2],spkT,'binWidth',boxSize/nBins,'smooth',smth);
                m2 = minions.rescaleData(map.z,eps,1);
                m2 = inpaint_nans(m2);
                
                m1all(:,:,cnt) = m1; % save rate map in a stack
                m2all(:,:,cnt) = m2;
                
                cnt = cnt+1;
            end
        end
    end
    
    m2all = minions.nan_pad3(m2all,15); % this prevents edge effects during rotations
    
    % rotation PV
    szmap = size(m1all,1);
    halfmap = round(szmap/2);
    
    rotCC = nan(szmap*2-1,szmap*2-1,(360/nDeg));
    shft_range = -1*(szmap-1):szmap-1;
    
    % rotate em
    for iRot = 1:(360/nDeg)
        %                     cc = nan(szmap-1); % center
        cc = nan(szmap*2-1); % full map
        m2_rot = imrotate3(m2all,(iRot-1)*nDeg,[0 0 1]);
        
        % chop extra bins created by rotating (the nan corners)
        szrot = size(m2_rot,1);
        cut = floor((szrot-szmap)/2);
        if mod(szrot,2)
            inds = cut+1:szrot-cut-1; % odd size remove extra one
        else
            inds = cut+1:szrot-cut;
        end
        m2_rot = m2_rot(inds,inds,:);
        
        % shift em
        for ix = shft_range
            for iy = shft_range
                m2_shift = circshift(m2_rot,[ix iy 0]);
                % wipe out pixels that wrapped around from other side
                if ix < 0
                    m2_shift(end-abs(ix):end,:,:) = nan;
                elseif ix > 0
                    m2_shift(1:ix,:,:) = nan;
                end
                if iy < 0
                    m2_shift(:,end-abs(iy):end,:) = nan;
                elseif iy > 0
                    m2_shift(:,1:iy,:) = nan;
                end
                if size(m1all) ~= size(m2_shift)
                    warning('a bad thing has happened')
                end
                % cut extra pixels (all nans) to save time
                remove = all(isnan(m2_shift),3);
                clear tmp1 tmp2
                for i = 1:size(m1all,3)
                    tmp = m1all(:,:,i);
                    tmp1(:,:,i) = tmp(~all(remove,2),~all(remove,1));
                    tmp = m2_shift(:,:,i);
                    tmp2(:,:,i) = tmp(~all(remove,2),~all(remove,1));
                end
                pv = analyses.populationVectorCorrelation(tmp1,tmp2,'rows','pairwise');
                pv(pv==0) = nan;
                npix = sum(~isnan(pv(:)));
                cc(ix+max(shft_range)+1,iy+max(shft_range)+1) = nanmean(pv(:)) .* npix; % this scales to number of pixels included
            end
        end
        
        rotCC(:,:,iRot) = cc ./ (size(m1all,1)*size(m1all,2)); % this completes the normalization
        fprintf('Rot %d deg ... ',(iRot-1)*nDeg)
    end
    fprintf('\n')
    
    % get max and min corr values for each rotation
    clear mn mx
    for i = 1:(360/nDeg)
        mx(i) = nanmax(rotCC(:,:,i),[],'all');
        mn(i) = nanmin(rotCC(:,:,i),[],'all');
    end
    
    % who wins?
    [~,mind] = nanmax(mx);
    fprintf('Rot = %d deg\n',(mind-1)*nDeg)
    
    % plot the ensemble CC at 0 degrees
    cc = rotCC(:,:,1);
    figure('InvertHardCopy', 'off');
    colorMapBRK(cc);
    colorbar
    hold on
    plot([szmap szmap],[0 szmap*2],'k-','linew',2)
    plot([0 szmap*2],[szmap szmap],'k-','linew',2)
    
%     saveas(gcf,fullfile(setFolder,sprintf('noRotationCC.jpg')))
%     saveas(gcf,fullfile(setFolder,sprintf('noRotationCC.fig')))
%     close(gcf)
    
    % plot the winner and note the central peak
    cc = rotCC(:,:,mind);
    figure('InvertHardCopy', 'off');
    colorMapBRK(cc);
    colorbar
    hold on
    plot([szmap szmap],[0 szmap*2],'k-','linew',2)
    plot([0 szmap*2],[szmap szmap],'k-','linew',2)
    
    try
        [r,c] = ind2sub(size(cc),find(imregionalmax(cc)));
    catch
        cc(isnan(cc)) = 0;
        [r,c] = ind2sub(size(cc),find(imregionalmax(cc)));
    end
    
    [ind,dist] = knnsearch([c,r],[size(cc,1)/2 size(cc,1)/2]);
    plot(c(ind),r(ind),'ko','markers',15,'linew',2)
    title(sprintf('Rotation: %d deg\nShift: %d, %d',(mind-1)*nDeg,c(ind)-szmap,(r(ind)-szmap)*-1))
    colorbar
    
%     saveas(gcf,fullfile(setFolder,sprintf('bestRotationCC.jpg')))
%     saveas(gcf,fullfile(setFolder,sprintf('bestRotationCC.fig')))
%     close(gcf)
    
    % plotting
    xVals = 0:(360/nDeg)-1;
    figure;
    plot(xVals,mx,'linewidth',1.5)
    set(gca,'xtickmode','manual','xtick',0:20:120,'xticklabels',0:60:360,'box','off','fontsize',10)
    ylim([0 inf])
    xlabel('Degrees')
    ylabel('PV correlation')
    
%     saveas(gcf,fullfile(setFolder,sprintf('rotationShift.jpg')))
%     saveas(gcf,fullfile(setFolder,sprintf('rotationShift.fig')))
%     close(gcf)
    
    % store values for rotation only
    mxStore(iComp,:,1) = mx;
    maxCorr(iComp,1) = nanmax(mx);
    minCorr(iComp,1) = nanmin(mx);
    rotAngle(iComp,1) = (mind-1)*3;
    mindStore(iComp,1) = mind;
    numStore(iComp,1) = size(units,2);
    rotCCstore(:,:,:,iComp) = rotCC;
    unitStore{iComp,1} = units;
end





