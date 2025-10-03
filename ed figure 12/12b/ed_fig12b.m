load rateMaps_sim1382
cellsToPlot = [103 132 177]; 
figure;
cnt = 1;
for iCell = 1:size(mapsA,3)
    mapA = general.smooth(mapsA(:,:,iCell),6.67);
    
    subplot(2,3,cnt)
    colorMapBRK(mapA);
    title(sprintf('Cell %d',cellsToPlot(iCell)))
    cnt = cnt + 1;
end
for iCell = 1:size(mapsB,3)
    mapB = general.smooth(mapsB(:,:,iCell),6.67);    
    subplot(2,3,cnt)
    colorMapBRK(mapB);
    title(sprintf('Cell %d',cellsToPlot(iCell)))
    cnt = cnt + 1;
end

%%
load rateMaps_sim4858
cellsToPlot = [120 176 77];
figure;
cnt = 1;
for iCell = 1:size(mapsA,3)
    mapA = general.smooth(mapsA(:,:,iCell),6.67);
    
    subplot(2,3,cnt)
    colorMapBRK(mapA);
    title(sprintf('Cell %d',cellsToPlot(iCell)))
    cnt = cnt + 1;
end
for iCell = 1:size(mapsB,3)
    mapB = general.smooth(mapsB(:,:,iCell),6.67);    
    subplot(2,3,cnt)
    colorMapBRK(mapB);
    title(sprintf('Cell %d',cellsToPlot(iCell)))
    cnt = cnt + 1;
end

%%
load rateMaps_sim7347

cellsToPlot = [132 71 193]; 
figure;
cnt = 1;
for iCell = 1:size(mapsA,3)
    mapA = general.smooth(mapsA(:,:,iCell),6.67);
    
    subplot(2,3,cnt)
    colorMapBRK(mapA);
    title(sprintf('Cell %d',cellsToPlot(iCell)))
    cnt = cnt + 1;
end
for iCell = 1:size(mapsB,3)
    mapB = general.smooth(mapsB(:,:,iCell),6.67);    
    subplot(2,3,cnt)
    colorMapBRK(mapB);
    title(sprintf('Cell %d',cellsToPlot(iCell)))
    cnt = cnt + 1;
end

%%
load rateMaps_sim3202

cellsToPlot = [16 105 192]; 
figure;
cnt = 1;
for iCell = 1:size(mapsA,3)
    mapA = general.smooth(mapsA(:,:,iCell),6.67);
    
    subplot(2,3,cnt)
    colorMapBRK(mapA);
    title(sprintf('Cell %d',cellsToPlot(iCell)))
    cnt = cnt + 1;
end
for iCell = 1:size(mapsB,3)
    mapB = general.smooth(mapsB(:,:,iCell),6.67);    
    subplot(2,3,cnt)
    colorMapBRK(mapB);
    title(sprintf('Cell %d',cellsToPlot(iCell)))
    cnt = cnt + 1;
end

%%
load rateMaps_sim4283

cellsToPlot = [124 185 190]; 
figure;
cnt = 1;
for iCell = 1:size(mapsA,3)
    mapA = general.smooth(mapsA(:,:,iCell),6.67);
    
    subplot(2,3,cnt)
    colorMapBRK(mapA);
    title(sprintf('Cell %d',cellsToPlot(iCell)))
    cnt = cnt + 1;
end
for iCell = 1:size(mapsB,3)
    mapB = general.smooth(mapsB(:,:,iCell),6.67);    
    subplot(2,3,cnt)
    colorMapBRK(mapB);
    title(sprintf('Cell %d',cellsToPlot(iCell)))
    cnt = cnt + 1;
end

%%
load rateMaps_sim390
cellsToPlot = [185 130 143]; 
figure;
cnt = 1;
for iCell = 1:size(mapsA,3)
    mapA = general.smooth(mapsA(:,:,iCell),6.67);
    
    subplot(2,3,cnt)
    colorMapBRK(mapA);
    title(sprintf('Cell %d',cellsToPlot(iCell)))
    cnt = cnt + 1;
end
for iCell = 1:size(mapsB,3)
    mapB = general.smooth(mapsB(:,:,iCell),6.67);    
    subplot(2,3,cnt)
    colorMapBRK(mapB);
    title(sprintf('Cell %d',cellsToPlot(iCell)))
    cnt = cnt + 1;
end

%%
load mxStore_rotSim

figure('position',[680,299,286,679]);
for iRow = 1:size(mxStoreAll,1)

    subplot(3,1,iRow)
    plot(1:120,mxStoreAll(iRow,:))
    ylim([0 0.9])
    [m,i] = nanmax(mxStoreAll(iRow,:));
    ylabel('PV correlation')
    xlabel('Rotation')
    set(gca,'xtick',0:20:120,'xticklabels',0:60:360)
end
