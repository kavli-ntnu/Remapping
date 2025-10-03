
% Standard analysis pipeline for Neuropixels recordings, one probe and one session at a time.
%
%   USAGE
%       analysis(folder,probeName,varargin)
%       folder              recording directory (e.g. N:\neuropixel\shared_data\26718_Coconut\2020-09-11_12-29-58)
%       probeName           probe name, e.g. 'MEC' or 'HP'
%                           these names must be stored in the animal's parent folder with the following format:
%                           \shared_data\26718_Coconut\_MEC_18408405672.sn
%                           \shared_data\26718_Coconut\_HP_18408403232.sn
%
%   OPTIONAL NAME-VALUE PAIRS
%       'session'           session number (default = 1)
%       'sessionMod'        one or more sessions that you want to call the session number given by 'session' (eg. [2] or [2, 3, 4], default = [])
%       'nBins'             number of spatial bins [x y] (default = [40 40], can give just 1 number)
%       'boxSize'           size of box [x y] in cm (or m) (default = [150 150], can give just 1 number)
%       'plot'              cell array or string specifying what to plot {'maps','tc','pps','pprainbow'}
%                           (default = all, skip plotting with 0, false, or, 'none')
%       'ksVer'             Kilosort version (string, default = '2.5')
%       'ksIDs'             1 will use KiloSort IDs for cluster numbers in figures, 0 will use array indices (default = 1)
%       'clrSpk'            color for spikes on path plot (default light theme = NTNU blue, default dark theme = yellow)
%       'clrTC'             color for HD tuning curve (default light theme = NTNU blue, default dark theme = yellow)
%       'skipSI'            1 will skip spatial information for speed (default = 0)
%       'skipGrid'          1 will skip grid stats for speed (default = 0)
%       'skipHD'            1 will skip HD stats for speed (default = 0). does not affect figures.
%       'phyGroups'         cell array or string specifying phy groups to be included {'good','mua','noise','unsorted'} (default = {'good','mua','unsorted'})
%       'customPath'        folder specifying custom path (e.g. 'Christy' will load/save from \26718_Coconut\2020-09-11_12-29-58\Christy\...) (default = '')
%       'force'             1 will delete previous results and redo it, 0 will load data and create just figures (default = 0).
%
%   DEPENDENCIES
%       BNT                 https://github.com/brkanter/BNT
%       hippo               https://github.com/brkanter/hippo
%       npx_pipeline        Contact Rich Gardner (richard[dot]gardner[at]ntnu[dot]no)
%
% Written by BRK 2020 (brkanter[at]gmail[dot]com)

function analysis(folder,probeName,varargin)

%% parse inputs
inp = inputParser();
inp.addRequired('folder')
inp.addRequired('probeName')
inp.addParameter('session',1,@(x) isnumeric(x) && numel(x) == 1);
inp.addParameter('sessionMod',[],@(x) isnumeric(x));
inp.addParameter('nBins',[40 40],@(x) isnumeric(x) && numel(x) <= 2);
inp.addParameter('boxSize',[150 150],@(x) isnumeric(x) && numel(x) <= 2);
inp.addParameter('plot',{'maps','tcs','pps','pprainbows'});
inp.addParameter('ksVer','2.5',@(x) ischar(x) || isstring(x) == 1);
% color theme depends on figure background color
ntnuBlue = [0,0.3086,0.6172];
clrSpk = ntnuBlue;
clrTC = ntnuBlue;
if sum(get(groot,'defaultfigurecolor')) < 0.7
    clrSpk = 'y';
    clrTC = 'y';
end
inp.addParameter('ksIDs',1,@(x) isnumeric(x) || islogical(x));
inp.addParameter('clrSpk',clrSpk);
inp.addParameter('clrTC',clrTC);
inp.addParameter('skipSI',0,@(x) isnumeric(x) || islogical(x));
inp.addParameter('skipGrid',0,@(x) isnumeric(x) || islogical(x));
inp.addParameter('skipHD',0,@(x) isnumeric(x) || islogical(x));
inp.addParameter('phyGroups',{'good','mua','unsorted'},@(x) all(ismember(lower(x),{'good','mua','noise','unsorted'})));
inp.addParameter('customPath','');
inp.addParameter('force',0,@(x) isnumeric(x) || islogical(x));
inp.parse(folder,probeName,varargin{:});
P = inp.Results;

% extract and curate
sessionNum = P.session;
sessionMod = P.sessionMod;
nBins = P.nBins;
if numel(nBins) == 1
    nBins = [nBins nBins];
end
boxSize = P.boxSize;
if numel(boxSize) == 1
    boxSize = [boxSize boxSize];
end
PLOT = P.plot;
if (isnumeric(PLOT) || islogical(PLOT)) && ~PLOT
    PLOT = 0;
end
phyGroupSelection = lower(P.phyGroups);
if ~iscell(phyGroupSelection)
    phyGroupSelection = {phyGroupSelection};
end
if isnumeric(probeName)
    probeName = num2str(probeName);
end
if ~isempty(P.customPath)
    customPath = ['\' P.customPath];
else
    customPath = P.customPath;
end

% check path exists
if isempty(dir(folder))
    error('Path does not exist. Check drive letter first.')
end

% match probe name to serial number
splits = regexp(folder,'\','split');
files = dir([strjoin(splits(1:end-1),'\'),'\*.sn']);
if numel(files)
    flag = true;
    i = 1;
    while flag
        if i > numel(files)
            error(sprintf('Failed to find .sn file for probe ''%s'' in the animal''s parent folder.',probeName))
        end
        if strfind(files(i).name,probeName)
            splits2 = strsplit(files(i).name,'_');
            if strcmpi(splits2{2},probeName)
                probeSN = splits2{end}(1:end-3);
                flag = false;
            else
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end
end

if any(boxSize < 40)
    warning(sprintf('You entered box size = %.2f x %.2f cm, but that can''t be right. Changing to %d x %d cm ...\n', ...
        boxSize(1),boxSize(2),boxSize(1)*100,boxSize(2)*100))
    boxSize = boxSize .* 100;
end

%% Load or run analysis
% If the data has been analyzed before, just (re)produce figures
if ~P.force && exist(sprintf('%s%s\\analysis_probe%s\\session%d\\results.mat',folder,customPath,probeName,sessionNum),'file')
    
    load(sprintf('%s%s\\analysis_probe%s\\session%d\\results.mat',folder,customPath,probeName,sessionNum),'results')
    fprintf('Loaded saved results from %s%s\\analysis_probe%s\\session%d\\results.mat ...\n',folder,customPath,probeName,sessionNum)
    
% Previous analysis not found, do everything
else

    fprintf('Processing data from %s ...\n',folder)
    fprintf('Probe ''%s'' (sn = %s) ...\n',probeName,probeSN)
    fprintf('Session %d ...\n',sessionNum)
    
    %% load data
  
    % note that we load the whole recording for now, extract some info, then narrow down to session
    % this is to ensure that results from each session are aligned
    d = npx_load(folder, "ksVersion", P.ksVer);
    if numel(d.sessions) < sessionNum
        error(sprintf('Session %d not found.',sessionNum))
    end
    
    % find probe with correct serial number
    flag = true;
    probeIdx = 1;
    while flag
        if probeIdx > numel(d.probes)
            error(sprintf('Hmmm...I didn''t find your probe ''%s'' :/',probeName))
        end
%         if strcmp(num2str(d.probes(probeIdx).metadata.ap.imProbeSN),probeSN)
        if strcmp(num2str(d.probes(probeIdx).metadata.ap.imDatPrb_sn),probeSN)
            flag = false;
        else
            probeIdx = probeIdx + 1;
        end
    end
    
    %% extract unit info
    N = numel(d.probes(probeIdx).units);
    ksID = [d.probes(probeIdx).units.id]';
    allSpikes = cell(N,1);
    phyGroup = cell(N,1);
    ksLabel = cell(N,1);
    for iUnit = 1:N
        allSpikes{iUnit} = d.probes(probeIdx).units(iUnit).spikeTimes;
        phyGroup{iUnit} = d.probes(probeIdx).units(iUnit).group;
        ksLabel{iUnit} = d.probes(probeIdx).units(iUnit).ks2Label;
    end
    
    %% sort by depth and filter by phy groups
    dpth = [d.probes(probeIdx).units.depth]';
    [dpthSrt,indSrt] = sort(dpth);
    allSpikes = allSpikes(indSrt);
    ksID = ksID(indSrt);
    phyGroup = phyGroup(indSrt);
    ksLabel = ksLabel(indSrt);
    keepDpth = ~isnan(dpthSrt); % depth is nan when there are 0 spikes so this removes dead clusters
    
    keepPhy = ismember(phyGroup,phyGroupSelection);
    keep = keepDpth & keepPhy;
    dpthSrt = dpthSrt(keep);
    allSpikes = allSpikes(keep);
    ksID = ksID(keep);
    phyGroup = phyGroup(keep);
    ksLabel = ksLabel(keep);
    N = numel(allSpikes);
    
    %% extract session
    buffer = 3; % sec, in case button press for session start/end was imprecise
    
    % check for renaming / merging sessions and also for restricted time ranges
    if ~isempty(sessionMod)
        tr = [];
        for iJoin = 1:length(sessionMod)
            tmp = d.sessions(sessionMod(iJoin)).validTimeRanges;
            tmp(1,1) = tmp(1,1) + buffer;
            tmp(end,2) = tmp(end,2) - buffer;
            tr = [tr; tmp];
        end
    else
        tr = d.sessions(sessionNum).validTimeRanges;
        tr(1,1) = tr(1,1) + buffer;
        tr(end,2) = tr(end,2) - buffer;
    end
    
    % only keep time range(s) of interest
    runInds = false(length(d.tracking.t),1);
    for iRange = 1:size(tr,1)
        [~, t0] = min(abs(d.tracking.t - tr(iRange,1)));
        [~, t1] = min(abs(d.tracking.t - tr(iRange,2)));
        runInds(t0:t1) = true;
    end

    %% pos and spikes only during task period
    posT = d.tracking.t(runInds);
    
    D = dir(folder);
    ind = find(~cellfun(@isempty,cellfun(@strfind,{D.name},repmat({'brokenTracking'},1,length(D)),'UniformOutput',false)));
    if isempty(ind)
        brokenTracking = false;
    else
        brokenTracking = true;
        posToriginal = posT;
        tok = strsplit(D(ind).name,'_');
        tok = str2double(tok{end}) / 1000;
        posT = posT - tok;
        fprintf('Found broken tracking file. Adjusted by %.3f secs ...\n',tok)
    end
    posX = d.tracking.x(runInds);
    posY = d.tracking.y(runInds);

    % make minor adjustments to tracking data to force consistent map sizes
    if abs((abs(range(posX)*100 - range(posY)*100)) - range(boxSize)) > 15
        % do not adjust more than 5 cm in either direction to avoid distortions
        warndlg(sprintf(['Tracking data does not match specified box size!\n\n' ...
            'X range: %.2f cm\nY range: %.2f cm\n' ...
            'Specified box size: %d x %d cm\n\n' ...
            'Skipping position adjustment so be aware map size could be unexpected.\n\n' ...
            '%s\nprobe%s\nsession%d'], ...
            range(posX)*100,range(posY)*100,boxSize(1),boxSize(2),folder,probeName,sessionNum))
        % just convert m to cm and start at (0,0)
        posX = posX + nanmin(posX);
        posX = posX * 100;
        posY = posY + nanmin(posY);
        posY = posY * 100;
    else
        fprintf('Rescaling position data to true box size %d x %d cm...\n',boxSize(1),boxSize(2))
        posX = minions.rescaleData(posX,0,boxSize(1));
        posY = minions.rescaleData(posY,0,boxSize(2));
    end
    
    % position data has been converted from meters to centimeters from this point forward
    
    HD = rad2deg(d.tracking.hd_azimuth(runInds));
    Fs = nanmedian(diff(posT));
    minTime = min(posT);
    maxTime = max(posT) + Fs;
    
    %% focus on active periods
    sm = 0; % best not to smooth position data
    x_sm = general.smooth(posX,sm);
    y_sm = general.smooth(posY,sm);
    t_sm = general.smooth(posT,sm);
    N_pos = size(posX,1);
    spd = zeros(N_pos,1);
    for i = 2:N_pos-1
        spd(i) = sqrt((x_sm(i+1) - x_sm(i-1))^2 + (y_sm(i+1) - y_sm(i-1))^2) / (t_sm(i+1) - t_sm(i-1));
    end
    spd(1) = spd(2);
    spd(end) = spd(end-1);
    smKer = 2*round(1/Fs); % 2 sec Gaussian kernel
    spd_sm = general.smooth(spd,smKer);
    
    % remove periods where running speed < 5 cm/s
    posX(spd_sm < 5) = nan;
    posY(spd_sm < 5) = nan;
    
    %% unit analysis
    spkPos = cell(N,1);
    spkInd = cell(N,1);
    maps = cell(N,1);
    meanRate = nan(N,1);
    peakRate = nan(N,1);
    si = nan(N,1);
    grids = nan(N,3);
    mvl = nan(N,1);
    meanAngle = nan(N,1);
    warning('off','BNT:numFields')
    
    c = 0; errGrid = 0; errHD = 0;
    for iUnit = 1:N
        
        c = c + 1;
        
        if iUnit == 1 || mod(iUnit,50) == 0
            fprintf('Unit %d of %d ...\n',c,N)
        end
        
        spkT = allSpikes{iUnit};
        spkT = spkT(spkT > minTime & spkT < maxTime);
        
        [spkPos{c},spkInd{c}] = data.getSpikePositions(spkT,[posT posX posY]);
        
        map = analyses.map([posT,posX,posY],spkT,'binWidth',boxSize./nBins,'smooth',2);
        maps(c) = {map.z};
        if ~isempty(spkPos{c})
            meanRate(c) = analyses.meanRate(spkPos{c}(:,1),[posT posX posY]);
        end
        if ~isempty(map.z)
            peakRate(c) = map.peakRate;
        end
        
        if ~P.skipSI
            infos = analyses.mapStatsPDF(map);
            si(c) = infos.content;
        end
        if ~P.skipGrid
            try
                acorr = analyses.autocorrelation(map.z);
                [gScore,gStats] = analyses.gridnessScore(acorr);
                % store as [score, spacing, orientation]
                grids(c,:) = [gScore, ...
                              nanmedian(gStats.spacing(:)) * (boxSize(1)/nBins(1)), ... % adjust spacing to get cm
                              gStats.orientation(1)]; 
            catch
                errGrid = errGrid + 1;
            end
        end
        if ~P.skipHD
            spkHD = HD(spkInd{c});
            try % maybe not enough spikes
                tc = analyses.turningCurve(spkHD,HD,Fs,'binWidth',6);
                tcStat = analyses.tcStatistics(tc,6,20);
                mvl(c) = tcStat.r;
                meanAngle(c) = tcStat.mean;
            catch
                errHD = errHD + 1;
            end
        end
        
    end
    warning('on','BNT:numFields')

    %% save results
    results.folder = folder;
    results.probeName = probeName;
    results.probeSN = probeSN;
    results.probeIdx = probeIdx;
    results.sessionNum = sessionNum;
    results.phyGroupSelection = phyGroupSelection;
    results.dpth = dpthSrt;
    results.ksID = ksID;
    results.phyGroup = phyGroup;
    results.ksLabel = ksLabel;
    results.N = N;
    results.posT = posT;
    results.posX = posX;
    results.posY = posY;
    results.maps = maps;
    results.meanRate = meanRate;
    results.peakRate = peakRate;
    results.spkPos = spkPos;
    results.spkInd = spkInd;
    results.hd = HD;
    if ~P.skipGrid
        results.grids = grids;
    end
    if ~P.skipSI
        results.si = si;
    end
    if ~P.skipHD
        results.mvl = mvl;
        results.meanAngle = meanAngle;
    end
    
    if ~exist(sprintf('%s%s\\analysis_probe%s',folder,customPath,probeName),'dir')
        mkdir(sprintf('%s%s\\analysis_probe%s',folder,customPath,probeName))
    end
    if ~exist(sprintf('%s%s\\analysis_probe%s\\session%d',folder,customPath,probeName,sessionNum),'dir')
        mkdir(sprintf('%s%s\\analysis_probe%s\\session%d',folder,customPath,probeName,sessionNum))
    end
    
    if errGrid
        warning('Grid stats failed for %d of %d units',errGrid,N)
    end
    if errHD
        warning('HD stats failed for %d of %d units',errHD,N)
    end
    fprintf('Saving results to %s%s\\analysis_probe%s\\session%d\\results.mat ...\n',folder,customPath,probeName,sessionNum)
    results = orderfields(results);
    save(sprintf('%s%s\\analysis_probe%s\\session%d\\results.mat',folder,customPath,probeName,sessionNum),'results')
    if ~isempty(sessionMod)
        save(sprintf('%s%s\\analysis_probe%s\\session%d\\sessionMod.mat',folder,customPath,probeName,sessionNum),'sessionMod')
    end
    
end

%% make figures
if ~isnumeric(PLOT)
    if find(strcmpi(PLOT,'maps'))
        fprintf('Making rate maps...\n')
        plotEm('maps',results,P.ksIDs,P.clrSpk,P.clrTC,customPath)
    end
    if find(strcmpi(PLOT,'tcs'))
        fprintf('Making tuning curves...\n')
        plotEm('tcs',results,P.ksIDs,P.clrSpk,P.clrTC,customPath)
    end
    if find(strcmpi(PLOT,'pps'))
        fprintf('Making path plots...\n')
        plotEm('pps',results,P.ksIDs,P.clrSpk,P.clrTC,customPath)
    end
    if find(strcmpi(PLOT,'pprainbows'))
        fprintf('Making HD path plot...\n')
        plotEm('pprainbows',results,P.ksIDs,P.clrSpk,P.clrTC,customPath)
    end
end
fprintf('Done!\n')

end


%% plotting function
function plotEm(type,results,ksIDs,clrSpk,clrTC,customPath)

pageWidth = 29.7; % A4 paper
pageHeight = 21;
spCols = 10;
spRows = 5;
leftEdge = 0;
rightEdge = 0;
topEdge = 0;
bottomEdge = 0;
spaceX = 0.2;
spaceY = 0;
sub_pos = subplot_pos(pageWidth,pageHeight,leftEdge,rightEdge,topEdge,bottomEdge,spCols,spRows,spaceX,spaceY);
figure('PaperUnits','cent','PaperSize',[pageWidth pageHeight],'PaperPos',[0 0 pageWidth pageHeight],'inverthardcopy','off','visible','off');
set(gcf,'color','w')
i = 1; pageNum = 1;

if strcmpi(type,'tcs')
    Fs = nanmedian(diff(results.posT));
end
if strcmpi(type,'pprainbows')
    hmap(1:360,1) = linspace(0,1,360);
    hmap(:,[2 3]) = 0.85; % brightness
    huemap = hsv2rgb(hmap);
    cmap = huemap;
end

% you want KiloSort IDs, but you cannot have them
if ksIDs && ~isfield(results,'ksID')
    warning('Did not find KS IDs, using cell indices instead.')
    ksIDs = false;
end

while i <= results.N
    
    if i == 1 || mod(i,50) == 0
        fprintf('Unit %d of %d ...\n',i,results.N)
    end
    for iRow = 1:spRows
        for iCol = 1:spCols
            
            axes('pos',sub_pos{iRow,iCol});
            
            if strcmpi(type,'maps')
                
                try
                    colorMapBRK(results.maps{i},'ydir','normal');
                    axis equal
                catch
                    axis off
                end
                if ksIDs
                    h = title(sprintf('%d',results.ksID(i)),'fontsize',14,'fontweight','bold');
                else
                    h = title(sprintf('%d',i),'fontsize',14,'fontweight','bold');
                end
                h.Position(2) = size(results.maps{i},1) + 2;
                
            elseif strcmpi(type,'tcs')
                
                if ksIDs
                    h = title(sprintf('%d',results.ksID(i)),'fontsize',14,'fontweight','bold');
                else
                    h = title(sprintf('%d',i),'fontsize',14,'fontweight','bold');
                end
                if numel(results.spkInd{i}) > 20
                    spkHD = results.hd(results.spkInd{i});
                    tc = analyses.turningCurve(spkHD,results.hd,Fs,'binWidth',6);
                    circularTurningBRK(tc(:,2)/max(tc(:,2)),'color',clrTC,'linewidth',3)
                    hold on
                    circularTurningBRK(tc(:,3)/max(tc(:,3)),'adjustaxis',false,'color',ones(1,3)*0.5)
                    hold off
                    if strcmpi(get(gca,'ydir'),'reverse')
                        set(gca,'ydir','normal')
                    end
                    h.Position(2) = 1.2;
                end
                h.Position(2) = 1.05;
                axis equal
                axis off
                
            elseif strcmpi(type,'pps')
                
                
                plot(results.posX,results.posY,'-','color',ones(1,3)*0.5)
                if numel(results.spkInd{i})
                    hold on
                    scatter(results.spkPos{i}(:,2),results.spkPos{i}(:,3),8,clrSpk,'filled','markerfacealpha',0.8);
                end
                axis equal
                axis off
                if ksIDs
                    h = title(sprintf('%d',results.ksID(i)),'fontsize',14,'fontweight','bold');
                else
                    h = title(sprintf('%d',i),'fontsize',14,'fontweight','bold');
                end
                h.Position(2) = max(results.posY(:)) + max(results.posY(:))*0.1;
                
            elseif strcmpi(type,'pprainbows')
                
                if numel(results.spkInd{i})
                    spkHD = results.hd(results.spkInd{i});
                    [vals,inds] = sort(spkHD);
                    spikePosSort = results.spkPos{i}(inds,:);
                    spkcmap = zeros(length(vals),3);
                    for iSpike = 1:length(vals)
                        try
                            spkcmap(iSpike,:) = cmap(round(vals(iSpike)),:);
                        catch
                            spkcmap(iSpike,:) = cmap(1,:);
                        end
                    end
                    
                    plot(results.posX,results.posY,'-','color',ones(1,3)*0.5)
                    hold on
                    scatter(spikePosSort(:,2),spikePosSort(:,3),8,spkcmap,'filled','markerfacealpha',0.8);
                    axis equal
                end
                axis off
                if ksIDs
                    h = title(sprintf('%d',results.ksID(i)),'fontsize',14,'fontweight','bold');
                else
                    h = title(sprintf('%d',i),'fontsize',14,'fontweight','bold');
                end
                h.Position(2) = max(results.posY(:)) + max(results.posY(:))*0.1;
                
            end
            
            if i == results.N
                set(gcf,'position',get(0,'screensize'));
                saveas(gcf,sprintf('%s%s\\analysis_probe%s\\session%d\\%s%d.png',results.folder,customPath,results.probeName,results.sessionNum,type,pageNum))
                close(gcf);
                return
            elseif mod(i,spCols*spRows) == 0
                set(gcf,'position',get(0,'screensize'));
                saveas(gcf,sprintf('%s%s\\analysis_probe%s\\session%d\\%s%d.png',results.folder,customPath,results.probeName,results.sessionNum,type,pageNum))
                close(gcf);
                figure('PaperUnits','cent','PaperSize',[pageWidth pageHeight],'PaperPos',[0 0 pageWidth pageHeight],'inverthardcopy','off','visible','off');
                pageNum = pageNum + 1;
            end
            
            i = i + 1;
            
        end
        
    end
    
end

end

