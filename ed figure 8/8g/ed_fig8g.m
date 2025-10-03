load pvStore_ed8g

figure('position',[-1566,54,1533,257]);
for iExp = 1:size(mxAA_CA1,1)
    subplot(1,4,2)
    plot(1:120,mxAA_CA1(iExp,:),'-','color',[1 153/255 153/255])
    hold on
    set(gca,'xticklabels',0:60:360)
    set(gca,'ytick',0:0.2:1)
    ylabel('PV correlation')
    title('CA1 AA')
    ylim([0 0.7])
end

for iExp = 1:size(mxAB_CA1,1)
    subplot(1,4,1)
    plot(1:120,mxAB_CA1(iExp,:),'-','color',[153/255 0 0])
    hold on
    set(gca,'xticklabels',0:60:360)
    set(gca,'ytick',0:0.2:1)
    ylabel('PV correlation')
    title('CA1 AB')
    ylim([0 0.7])
end

for iExp = 1:size(mxAA_CA3,1)
    subplot(1,4,4)
    plot(1:120,mxAA_CA3(iExp,:),'-','color',[153/255 204/255 1])
    hold on
    set(gca,'xticklabels',0:60:360)
    set(gca,'ytick',0:0.2:1)
    ylabel('PV correlation')
    title('CA3 AA')
    ylim([0 0.7])
end

for iExp = 1:size(mxAB_CA3,1)
    subplot(1,4,3)
    plot(1:120,mxAB_CA3(iExp,:),'-','color',[0 102/255 204/255])
    hold on
    set(gca,'xticklabels',0:60:360)
    set(gca,'ytick',0:0.2:1)
    ylabel('PV correlation')
    title('CA3 AB')
    ylim([0 0.7])
end

