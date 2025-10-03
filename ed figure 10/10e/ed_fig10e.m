load fig3de_sessionPairs
load expList_HP_3de
groupIDs = [1 3 2];

load colorMat
load mxStore_ed10e

figure('position',[-1714,23,1614,282]);
for iGrp = 1:size(groupIDs,2)
    subplot(1,4,iGrp)
    inds = find(idx == groupIDs(iGrp));
    mxTemp = mxStoreAll(inds,:);

    for iRow = 1:size(mxTemp,1)
        if inds(iRow) == 23
            m = mxTemp(iRow,:);
            [maxVal,ind] = nanmax(m);
            m = circshift(m,-ind);
            plot(1:120,m,'color',colorMat(iGrp,:))
        else
        
        plot(1:120,mxTemp(iRow,:),'color',colorMat(iGrp,:))
        hold on
        end
    end
    ylim([0 0.75])
    set(gca,'ytick',0:0.25:1,'xtick',0:20:120,'xticklabels',0:60:360)
    ylabel('PV correlation')
    xlabel('Rotation')
end

load mxStable_ed10e
subplot(1,4,4)

for iRow = 1:size(mxStable,1)   
    plot(1:120,mxStable(iRow,:),'color',[0.5 0.5 0.5])
    hold on
end
ylim([0 0.75])
set(gca,'ytick',0:0.25:1,'xtick',0:20:120,'xticklabels',0:60:360)
ylabel('PV correlation')
xlabel('Rotation')
    
