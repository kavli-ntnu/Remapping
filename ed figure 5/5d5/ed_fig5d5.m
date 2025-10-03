
load distStore_M2_ed5d

figure('position',[682,305,370,314]);
for iRow = 1:size(distStore,1)
    plot(1:3,distStore(iRow,:),'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:3,distStore(iRow,:),'o','color',[0.5 0.5 0.5],'linewidth',1.5)
end

averageDist = nanmean(distStore,2);
plot(repmat(4,[size(averageDist,1) 1]),averageDist,'o','color',[0.5 0.5 0.5],'linewidth',1.5)

plot(1:3,nanmean(distStore(:,1:3)),'r+','markersize',15,'linewidth',2)
plot(4,nanmean(averageDist),'r+','markersize',15,'linewidth',2)

xlim([0.5 4.5])
ylim([0 90])
set(gca,'xtick',[1:4],'xticklabels',{'M1-M2','M1-M3','M2-M3','Mean'})
xtickangle(45)
set(gca,'ytick',0:20:120)
ylabel('Distance (cm)')
box on

%%
load randDist_rotM2
randDistMean = nan(size(finalDistStore,1),1);
for iExp = 1:size(randDistMean,1)
    randDistMean(iExp,1) = nanmean(nanmean(finalDistStore{iExp,1},2));
end
randDist = randDistMean;

plot(repmat(4.5,[size(randDist,1) 1]),randDist,'o','color',[0.5 0.5 0.5],'linewidth',1.5)
xlim([0.5 5])
plot(4.5,nanmean(randDist),'r+','markersize',15,'linewidth',2)
set(gca,'xtick',[1:4 4.5],'xticklabels',{'M1-M2','M1-M3','M2-M3','Mean','Shuffle'})

[h,p,stats] = ttest2(randDist,averageDist)

