%%
load distStoreM4

figure('position',[682,305,370,314]);
plot(1:6,distStore(1,:),'-','color',[0.8 0.8 0.8],'linewidth',1.5)
hold on
plot(1:6,distStore(1,:),'o','color',[0.5 0.5 0.5],'linewidth',1.5)

averageDist = nanmean(distStore,2);
plot(repmat(7,[size(averageDist,1) 1]),averageDist,'o','color',[0.5 0.5 0.5],'linewidth',1.5)

xlim([0.5 7.5])
ylim([0 90])
set(gca,'xtick',[1:7],'xticklabels',{'M1-M2','M1-M3','M1-M4','M2-M3','M2-M4','M3-M4','Mean'})
xtickangle(45)
set(gca,'ytick',0:20:120)
ylabel('Distance (cm)')
box on


