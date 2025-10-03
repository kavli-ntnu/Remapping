load borderMaps_ed11e

counter = 1;
page = 1;
numPlots = 4; 
figure('position',[1,41,1920,963],'inverthardcopy','off')
numRows = 6;

for iCell = 1:size(mapStore,1) 
    for iSession = 1:4
        if ~isempty(mapStore{iCell,iSession})
            map = mapStore{iCell,iSession};

            subplot(numRows,numPlots,counter)
            colorMapBRK(mapStore{iCell,iSession});
            text(0,45,sprintf('%.2f Hz',nanmax(nanmax(map))),'fontsize',8)
        end
        counter = counter + 1;
    end
end

