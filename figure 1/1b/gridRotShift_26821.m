
load expList_26821
load sData_MEC_26821
load gridIDstore_26821
numMods = size(gridIDstore,2);

iExp = 1;
sessionComp = [1 2];

%% set box size and nDeg
boxSize = expList_MEC{iExp,7}(1);
nBins = 40;
nDeg = 3;

%%
% initialize storage variables for total number of modules
mxStore = nan(numMods,360/nDeg);
maxCorr = nan(1,numMods);
rotAngle = nan(1,numMods);
mindStore = nan(1,numMods);
shiftStore = nan(numMods,2);
numStore = nan(1,numMods);
rotCCstore = nan(nBins*2-1,nBins*2-1,360/nDeg,numMods);
unitStore = cell(1);

for iMod = 1:numMods
    
    units_probe1 = gridIDstore{1,iMod}';
    units_probe2 = gridIDstore{2,iMod}';
    units = units_probe1;
    
    %% make map stacks
    smth = 1;
    
    % rescale pos data to fit box and to match each other
    posT1 = sData{1,sessionComp(1)}.results.posT; % timestamps for X and Y
    posX1 = minions.rescaleData(sData{1,sessionComp(1)}.results.posX,0,boxSize); % positions X
    posY1 = minions.rescaleData(sData{1,sessionComp(1)}.results.posY,0,boxSize); % positions Y
    min1x = nanmin(posX1);
    max1x = nanmax(posX1);
    min1y = nanmin(posY1);
    max1y = nanmax(posY1);
    posT2 = sData{1,sessionComp(2)}.results.posT;
    posX2 = minions.rescaleData(sData{1,sessionComp(2)}.results.posX,min1x,max1x);
    posY2 = minions.rescaleData(sData{1,sessionComp(2)}.results.posY,min1y,max1y);
    
    cnt = 1;
    clear m1all m2all
    
    % maps are rescaled from 0 to 1 to ignore rate diffs and nans are filled in
    for iUnit = units
        
        if ~isempty(sData{1,sessionComp(1)}.results.spkPos{iUnit})
            spkT = sData{1,sessionComp(1)}.results.spkPos{iUnit}(:,1); % get spike positions for this cell in session 1
            map = analyses.map([posT1,posX1,posY1],spkT,'binWidth',boxSize/nBins,'smooth',smth); % plot position data with spike times
            m1 = minions.rescaleData(map.z,eps,1); % m1 = rescaled rate map for this cell in first session, avoid having 0s in map
            m1 = inpaint_nans(m1); % interpolate surrounding values to fill in nans
            
            if ~isempty(sData{1,sessionComp(2)}.results.spkPos{iUnit})
                spkT = sData{1,sessionComp(2)}.results.spkPos{iUnit}(:,1);
                map = analyses.map([posT2,posX2,posY2],spkT,'binWidth',boxSize/nBins,'smooth',smth);
                m2 = minions.rescaleData(map.z,eps,1);
                m2 = inpaint_nans(m2);
                
                m1all(:,:,cnt) = m1; % save rate map in a stack
                m2all(:,:,cnt) = m2;
                
                cnt = cnt+1;
            end
        end
    end
    
    units = units_probe2;
    for iUnit = units
        if ~isempty(sData{2,sessionComp(1)}.results.spkPos{iUnit})
            spkT = sData{2,sessionComp(1)}.results.spkPos{iUnit}(:,1); % get spike positions for this cell in session 1
            map = analyses.map([posT1,posX1,posY1],spkT,'binWidth',boxSize/nBins,'smooth',smth); % plot position data with spike times
            m1 = minions.rescaleData(map.z,eps,1); % m1 = rescaled rate map for this cell in first session, avoid having 0s in map
            m1 = inpaint_nans(m1); % interpolate surrounding values to fill in nans
            
            if ~isempty(sData{2,sessionComp(2)}.results.spkPos{iUnit})
                spkT = sData{2,sessionComp(2)}.results.spkPos{iUnit}(:,1);
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
    
    %% module xcorrs
    
    szmap = size(m1all,1);
    halfmap = round(szmap/2);
    
    % only shift halfway (center)
    %         shft_range = -1*(halfmap-1):halfmap-1;
    %         rotCC = nan(szmap-1,szmap-1,(360/nDeg));
    
    % shift all the way to the edge (full map)
    shft_range = -1*(szmap-1):szmap-1;
    rotCC = nan(szmap*2-1,szmap*2-1,(360/nDeg));
    
    % rotate em
    for iRot = 1:(360/nDeg)
        cc = nan(szmap-1); % center
        %     cc = nan(szmap*2-1); % full map
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
                % OK to turn nans to zeros here because nans are matched for all cells in stack (i.e. the whole pixel), so corr
                % will yield nan
                tmp1(isnan(tmp1)) = 0;
                tmp2(isnan(tmp2)) = 0;
                pv = analyses.populationVectorCorrelation(tmp1,tmp2,'rows','complete');
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
    
    % show each rotation with same color scale
    figure('InvertHardCopy', 'off');
    set(gcf,'color','w')
    for i = 1:(360/nDeg)
        subplot(ceil(sqrt(360/nDeg)),ceil(sqrt(360/nDeg)),i)
        colorMapBRK(rotCC(:,:,i),'cutoffs',[min(mn) max(mx)]);
        title(sprintf('%.2f',nanmax(rotCC(:,:,i),[],'all')))
    end
%     close(gcf)
    
    % who wins?
    [~,mind] = nanmax(mx);
    fprintf('Rot = %d deg\n',(mind-1)*nDeg)
    
    % plot the winner and note the central peak
    cc = rotCC(:,:,mind);
    figure('InvertHardCopy', 'off');
    set(gcf,'color','w')
    colorMapBRK(cc);
    hold on
    %         plot([szmap/2 szmap/2],[0 szmap],'k-','linew',2)
    %         plot([0 szmap],[szmap/2 szmap/2],'k-','linew',2)
    plot([szmap szmap],[0 szmap*2],'k-','linew',2)
    plot([0 szmap*2],[szmap szmap],'k-','linew',2)
        
%     saveas(gcf,fullfile(outFolder,'bestModRotation.jpg'))
%     close(gcf)
%     fprintf('Pics saved here: %s\n',outFolder)
    
    %% plot correlation by angle
    xVals = 0:nDeg:359;
    figure;
    set(gcf,'color','w')
    plot(xVals,mx)
    set(gca,'xtickmode','manual','xtick',0:30:360,'xticklabels',0:60:360,'box','off','fontsize',10,'fontweight','bold')
    % ylim([0 inf])
    xlabel('Degrees')
    ylabel('Correlation')
%     saveas(gcf,fullfile(outFolder,'correlationDeg.jpg'))
%     close(gcf)
    
    %% store values for this module before moving on
    mxStore(iMod,:) = mx;
    maxCorr(1,iMod) = nanmax(mx);
    rotAngle(1,iMod) = (mind-1)*3;
    mindStore(1,iMod) = mind;
    shiftStore(iMod,:) = [c(ind)-szmap/2 r(ind)-szmap/2];
    numStore(1,iMod) = length(units);
    rotCCstore(:,:,:,iMod) = rotCC;
    unitStore{1,iMod} = units;
end

% save summary results of all modules
% save(fullfile(outFolder,'gridRotShift.mat'),'mxStore','maxCorr','rotAngle','mindStore','shiftStore','numStore','rotCCstore','unitStore','boxSize','nBins','nDeg')

%% plot correlation values for each module in the same plot
xVals = 0:nDeg:359;
figure;
set(gcf,'color','w')
hold on

for iMod = 1:numMods
    plot(xVals,mxStore(iMod,:))
    plot(rotAngle(1,iMod),maxCorr(1,iMod),'k*')
    %         xline(rotAngle(1,iMod))
end
% xlim([-1 360])

set(gca,'xtickmode','manual','xtick',0:30:360,'xticklabels',0:30:360,'box','off','fontsize',10)
% ylim([0 inf])
xlabel('Degrees')
ylabel('Correlation')

% saveas(gcf,fullfile(outFolder,'rotationByMod.jpg'))
% close(gcf)

