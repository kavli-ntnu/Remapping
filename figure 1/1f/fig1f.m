load expList_MEC_1f
load pkStoreCompiled_1f

load colorMapLarge
colorMap(10,:) = [0 0 0];

figure('position',[-1944,327,1717,142]);
cnt = 1;
for iCol = 1:size(pkStoreCompiled,2)
    temp = pkStoreCompiled{1,iCol};

    subplot(1,11,cnt)
    for iMod = 1:size(temp,1)
        a = plot([40 temp(iMod,1)],[40 temp(iMod,2)],'-','linewidth',2,'color',colorMap(iMod,:));
        hold on
%         plot(temp(iMod,1),temp(iMod,2),'.','color',colorMap(iMod,:),'markersize',15)
        set(gca,'ydir','rev')
        xlim([20 60])
        ylim([20 60])
        axis square
        set(gca,'xtick','','ytick','')
        title(sprintf('%s',expList_MEC{iCol,2}(1,1:5)))
        xline(40)
        yline(40)
    end
    cnt = cnt + 1;
end
