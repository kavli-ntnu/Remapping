load expList_26718
load rotCCstore_26718

numMods = size(rotCCstore,4);
     
%%
figure;
cnt = 1;
for iMod = 1:3

    szmap = 40;
    indices = [21 31 38 39 40 41 42 43 44 51 61]; 

    maxCorr = mxStore(iMod,mindStore(iMod));
    
    for ind = 1:size(indices,2)
        subplot(3,size(indices,2),cnt)
        colorMapBRK(rotCCstore(:,:,indices(ind),iMod),'cutoffs',[-inf maxCorr]);
        hold on
        plot([szmap szmap],[0 szmap*2],'k-','linew',1)
        plot([0 szmap*2],[szmap szmap],'k-','linew',1)
    
        title(sprintf('%.2f',mxStore(iMod,indices(ind))))
        cnt = cnt + 1;
    end
end
      
%%
load colorMap

figure;
for iMod = 1:numMods
    plot(1:120,mxStore(iMod,:),'color',colorMap(iMod,:),'linewidth',1.5)
    hold on
end
box on
xlabel('Rotation (deg)')
ylabel('Correlation')
set(gca,'xtick',0:20:120,'xticklabels',0:60:360)
set(gca,'ytick',0:0.1:1)
