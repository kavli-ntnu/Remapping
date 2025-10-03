
load gridMaps_ed11d

counter = 1;
page = 1;
numPlots = 5; 
figure('position',[1,41,1920,963],'inverthardcopy','off')
numRows = 6;

for iCell = 1:size(mapStore,1) 
    for iSession = 1:4
        if ~isempty(mapStore{iCell,iSession})
            map = mapStore{iCell,iSession};
            
            acorr = analyses.autocorrelation(map);
            [gScore,gStats,~,COM] = analyses.gridnessScore(acorr);
            
            subplot(numRows,numPlots,counter)
            colorMapBRK(mapStore{iCell,iSession},'ydir','normal');
            text(0,45,sprintf('%.2f Hz',nanmax(nanmax(map))),'fontsize',8)
            
            if ~isempty(gScore)
                text(30,45,sprintf('%.2f',gScore),'fontsize',8)
            end
        end
        counter = counter + 1;
    end

    % to plot cc, data must be loaded in from singleCells.mat
    subplot(numRows,numPlots,counter)
    colorMapBRK(ccStore(:,:,iCell));
    hold on
    xline(39)
    yline(39)
    text(0,-10,sprintf('%.1f°',rotCell(iCell,1)),'fontsize',8)    
    counter = counter + 1;
end

