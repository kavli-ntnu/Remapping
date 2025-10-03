
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
simplePVstoreRaw = nan(nBins*nBins,numComps,1);

for iComp = 1:numComps
    
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
    clear m1all m2all pv
    
    % maps are rescaled from 0 to 1 to ignore rate diffs and nans are filled in
    for iUnit = units
        
        if ~isempty(sData{1,sessionComp(iComp,1)}.results.spkPos{iUnit})
            spkT = sData{1,sessionComp(iComp,1)}.results.spkPos{iUnit}(:,1); % get spike positions for this cell in session 1
            map = analyses.map([posT1,posX1,posY1],spkT,'binWidth',boxSize/nBins,'smooth',smth); % plot position data with spike times
            m1 = map.z;
       
            
            if ~isempty(sData{1,sessionComp(iComp,2)}.results.spkPos{iUnit})
                spkT = sData{1,sessionComp(iComp,2)}.results.spkPos{iUnit}(:,1);
                map = analyses.map([posT2,posX2,posY2],spkT,'binWidth',boxSize/nBins,'smooth',smth);  
                m2 = map.z;
  
                m1all(:,:,cnt) = m1; % save rate map in a stack
                m2all(:,:,cnt) = m2;
                
                cnt = cnt+1;
            end
        end
    end
    
    % simple PV (pixels)
    pv = analyses.populationVectorCorrelation(m1all,m2all,'rows','pairwise');
    simplePVstoreRaw(:,iComp,1) = reshape(pv,[],1);    
end

%%
figure;
hold on
colors = {'b','r','k','m'};  % A1-A2, B1-B2, A1-B1, B2-A2
for iComp = 1:numComps
    h = cdfplot(simplePVstoreRaw(:,iComp,1)); 
    set(h,'color',colors{iComp},'linewidth',2)
end
set(gca,'xtickmode','manual','box','off','fontsize',12)
set(gca,'ytick',0:0.2:1)
xlabel('PV correlation','fontsize',14,'fontweight','bold')
ylabel('Proportion','fontsize',14,'fontweight','bold')
title(sprintf('n = %d',size(units,2)))

nanmedian(simplePVstoreRaw)