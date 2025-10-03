load minDistStore_ed5d
distStore = minDistStore;

load minDistRand_ed5d

%%
figure('position',[682,305,370,314]);
for iRow = 1:size(distStore,1)
    thisRow = distStore(iRow,:);

    plot(1:size(thisRow,2),thisRow,'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:size(thisRow,2),thisRow,'o','color',[0.5 0.5 0.5],'linewidth',1.5)  
end

averageDist = nanmean(distStore,2);
plot(repmat(4,[size(averageDist,1) 1]),averageDist,'o','color',[0.5 0.5 0.5],'linewidth',1.5)

plot(1:3,nanmean(distStore(:,1:3)),'r+','markersize',15,'linewidth',2)
plot(4,nanmean(averageDist),'r+','markersize',15,'linewidth',2)
    
rDist = randDistExp(~isnan(randDistExp));
plot(repmat(4.5,[size(rDist,1) 1]),rDist,'o','color',[0.5 0.5 0.5],'linewidth',1.5)
plot(4.5,nanmean(rDist),'r+','markersize',15,'linewidth',2)
ylim([-1 45])

xlim([0.5 5])
set(gca,'xtick',[1:4 4.5],'xticklabels',{'M1-M2','M1-M3','M2-M3','Mean','Shuffle'})
xtickangle(45)
set(gca,'ytick',0:10:100)
ylabel('Distance (cm)')
box on

%%
load pvStore_ed5d
distStore = pvStore;

load pvStoreRand_ed5d

%%
figure('position',[682,305,370,314]);
for iRow = 1:size(distStore,1)
    thisRow = distStore(iRow,:);

    plot(1:size(thisRow,2),thisRow,'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:size(thisRow,2),thisRow,'o','color',[0.5 0.5 0.5],'linewidth',1.5)  
end

averageDist = nanmean(distStore,2);
plot(repmat(4,[size(averageDist,1) 1]),averageDist,'o','color',[0.5 0.5 0.5],'linewidth',1.5)

plot(1:3,nanmean(distStore(:,1:3)),'r+','markersize',15,'linewidth',2)
plot(4,nanmean(averageDist),'r+','markersize',15,'linewidth',2)
    
rDist = randDistPV;
plot(repmat(4.5,[size(rDist,1) 1]),rDist,'o','color',[0.5 0.5 0.5],'linewidth',1.5)
plot(4.5,nanmean(rDist),'r+','markersize',15,'linewidth',2)
ylim([-1 55])

xlim([0.5 5])
set(gca,'xtick',[1:4 4.5],'xticklabels',{'M1-M2','M1-M3','M2-M3','Mean','Shuffle'})
xtickangle(45)
set(gca,'ytick',0:10:100)
ylabel('Distance (cm)')
box on
