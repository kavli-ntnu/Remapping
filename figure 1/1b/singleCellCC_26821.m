
load expList_26821
load gridIDstore_26821
numMods = size(gridIDstore,2);
load sData_MEC_26821

sessionComp = [1 2];
iExp = 1;
sessionNums = expList_MEC{iExp,5};

%% set box size and nDeg
boxSize = expList_MEC{iExp,7}(1);
nBins = 40;
nDeg = 3;

%% set field detection thresholds
fieldThresh = 0.3;
binWidth = 5;
minBins = 8;
minPeak = 1;
threshRate = 7;

%%
for iProbe = 1 %:2
    
    maps = cell(sData{iProbe,2}.results.N,length(sessionNums));
    for iSession = sessionNums
        maps(:,iSession) = sData{iProbe,iSession}.results.maps;
    end
    
    %%
    for iMod = 1:numMods
        
        units = gridIDstore{iProbe,iMod}';
        
        %% make map stacks
        smth = 1;
        
        % rescale pos data to fit box and to match each other
        posT1 = sData{iProbe,sessionComp(1)}.results.posT; % timestamps for X and Y
        posX1 = minions.rescaleData(sData{iProbe,sessionComp(1)}.results.posX,0,boxSize); % positions X
        posY1 = minions.rescaleData(sData{iProbe,sessionComp(1)}.results.posY,0,boxSize); % positions Y
        min1x = nanmin(posX1);
        max1x = nanmax(posX1);
        min1y = nanmin(posY1);
        max1y = nanmax(posY1);
        posT2 = sData{iProbe,sessionComp(2)}.results.posT;
        posX2 = minions.rescaleData(sData{iProbe,sessionComp(2)}.results.posX,min1x,max1x);
        posY2 = minions.rescaleData(sData{iProbe,sessionComp(2)}.results.posY,min1y,max1y);
        
        cnt = 1;
        clear m1all m2all
        
        % maps are rescaled from 0 to 1 to ignore rate diffs and nans are filled in
        for iUnit = units
            
            if ~isempty(sData{iProbe,sessionComp(1)}.results.spkPos{iUnit})
                spkT = sData{iProbe,sessionComp(1)}.results.spkPos{iUnit}(:,1); % get spike positions for this cell in session 1
                map = analyses.map([posT1,posX1,posY1],spkT,'binWidth',boxSize/nBins,'smooth',smth); % plot position data with spike times
                m1 = minions.rescaleData(map.z,eps,1); % m1 = rescaled rate map for this cell in first session, avoid having 0s in map
                m1 = inpaint_nans(m1); % interpolate surrounding values to fill in nans
                
                if ~isempty(sData{iProbe,sessionComp(2)}.results.spkPos{iUnit})
                    spkT = sData{iProbe,sessionComp(2)}.results.spkPos{iUnit}(:,1);
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
        
        %% single cell xcorrs
        
        szmap = size(m1all,1);
        halfmap = round(szmap/2);
        
        % only shift halfway (center)
        % shft_range = -1*(halfmap-1):halfmap-1;
        
        % shift all the way to the edge (full map)
        shft_range = -1*(szmap-1):szmap-1;
        
        mindMod = mindStore(iMod); % loads module mind just obtained above
        
        numCells = numel(units);
        
        resCell = nan(numCells,10);
        mxStoreCell = nan(numCells,360/nDeg);
        
        ccStoreBest = nan(szmap*2-1,szmap*2-1,numCells);
        ccStoreAll = nan(szmap*2-1,szmap*2-1,360/nDeg,numCells);
        
        for iUnit = 1:numCells
            
            % get maps
            m1 = m1all(:,:,iUnit);
            m2 = m2all(:,:,iUnit);
            m2 = minions.nan_pad2(m2,15); % prevent edge effects during rotations
            
            % xcorr with best rotation for module
            % cc = nan(szmap-1); % center
            cc = nan(szmap*2-1); % full map
            
            m2_rot = imrotate(m2,(mindMod-1)*nDeg,'bilinear');
            szrot = size(m2_rot,1);
            cut = floor((szrot-szmap)/2);
            if mod(szrot,2)
                inds = cut+1:szrot-cut-1; % odd size remove extra one
            else
                inds = cut+1:szrot-cut;
            end
            m2_rot = m2_rot(inds,inds,:);
            
            % this looks like the module one above, but it's in 2D instead of 3D
            for ix = shft_range
                for iy = shft_range
                    m2_shift = circshift(m2_rot,[ix iy 0]);
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
                    if size(m1) ~= size(m2_shift)
                        warning('a bad thing has happened')
                    end
                    npix = min([sum(~isnan(m1(:))),sum(~isnan(m2_shift(:)))]);
                    
                    % use spatialCorr for 2D instead of pvCorr for 3D
                    cc(ix+max(shft_range)+1,iy+max(shft_range)+1) = analyses.spatialCrossCorrelation(m1,m2_shift) * npix;
                end
            end
            cc = cc ./ sum(~isnan(m1(:))); % normalization
            
            cc(isnan(cc)) = 0;
            [r,c] = ind2sub(size(cc),find(imregionalmax(cc)));
            [ind,dist] = knnsearch([c,r],[size(cc,1)/2 size(cc,1)/2]);
            
            resCell(iUnit,1) = units(iUnit); % unit #
            resCell(iUnit,2) = (mindMod-1)*nDeg; % best mod rotation
            resCell(iUnit,3) = c(ind); % x shift at mod best
            resCell(iUnit,4) = r(ind); % y shift at mod best
            
            ccStoreBest(:,:,iUnit) = cc;
            
            %%% xcorr with best rotation for unit %%%
            
            % rotCC = nan(szmap-1,szmap-1,(360/nDeg)); % center
            rotCC = nan(szmap*2-1,szmap*2-1,(360/nDeg)); % full map
            
            for iRot = 1:(360/nDeg)
                % cc = nan(szmap-1); % center
                cc = nan(szmap*2-1); % full map
                
                m2_rot = imrotate(m2,(iRot-1)*nDeg,'bilinear');
                szrot = size(m2_rot,1);
                cut = floor((szrot-szmap)/2);
                if mod(szrot,2)
                    inds = cut+1:szrot-cut-1; % odd size remove extra one
                else
                    inds = cut+1:szrot-cut;
                end
                m2_rot = m2_rot(inds,inds,:);
                for ix = shft_range
                    for iy = shft_range
                        m2_shift = circshift(m2_rot,[ix iy 0]);
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
                        if size(m1) ~= size(m2_shift)
                            warning('a bad thing has happened')
                        end
                        npix = min([sum(~isnan(m1(:))),sum(~isnan(m2_shift(:)))]);
                        cc(ix+max(shft_range)+1,iy+max(shft_range)+1) = analyses.spatialCrossCorrelation(m1,m2_shift) * npix;
                    end
                end
                rotCC(:,:,iRot) = cc ./ sum(~isnan(m1(:)));
            end
            
            for i = 1:(360/nDeg)
                mx(i) = nanmax(rotCC(:,:,i),[],'all');
                mn(i) = nanmin(rotCC(:,:,i),[],'all');
            end
            [~,mind] = nanmax(mx);
            mxStoreCell(iUnit,:) = mx;
            
            cc(isnan(cc)) = 0;
            [r,c] = ind2sub(size(cc),find(imregionalmax(cc)));
            [ind,dist] = knnsearch([c,r],[size(cc,1)/2 size(cc,1)/2]);
            
            resCell(iUnit,5) = (mind-1)*nDeg; % best cell rotation
            resCell(iUnit,6) = c(ind); % x shift
            resCell(iUnit,7) = r(ind); % y shift
            
            ccStoreAll(:,:,:,iUnit) = rotCC;
        end       
    end
end

