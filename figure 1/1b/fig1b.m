load expList_26821
load gridCellCC_26821
load sData_MEC_26821

gridCells = [523 538 373 388 147 224]; % grid cells to plot
sessionNums = expList_MEC{1,5};

%% extract maps
maps = cell(sData{1,1}.results.N,size(sessionNums,2));
si = nan(sData{1,1}.results.N,size(sessionNums,2));
depth = nan(sData{1,1}.results.N,size(sessionNums,2));

for iSession = sessionNums
    maps(:,iSession) = sData{1,iSession}.results.maps;
    si(:,iSession) = sData{1,iSession}.results.si;
    depth(:,iSession) = sData{1,iSession}.results.dpth;
end

%%
numRows = 6;
counter = 1;
page = 1;
numPlots = 3; % (2 rate maps, 1 CC)
newPg = 1:6*numPlots:size(gridCells,1)*3;

for iCell = 1:size(gridCells,2)
    if sum(counter == newPg)
        figure('position',[1,41,1920,963],'inverthardcopy','off')
        counter = 1;
    end
    
    for iSession = 1:2
        if ~isempty(maps{gridCells(iCell),iSession})
            map = maps{gridCells(iCell),iSession};
          
            subplot(numRows,numPlots,counter)
            colorMapBRK(maps{gridCells(iCell),iSession},'ydir','normal');
            text(0,65,sprintf('%d Hz',round(nanmax(nanmax(map)))),'fontsize',8)
        end
        counter = counter + 1;
        
        if iSession == 1
            text(-25,20,sprintf('#%d',gridCells(iCell)),'fontsize',8)
        end
    end
    
    subplot(numRows,numPlots,counter)
    colorMapBRK(ccStore(:,:,iCell));
    hold on
    xline(size(ccStore,1)/2)
    yline(size(ccStore,1)/2)
    
    counter = counter + 1;
end