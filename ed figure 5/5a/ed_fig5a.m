load expList_MEC_1f

load normDistAA_ed5a
normDist = normDistAA;

figure('position',[723,214,274,314]);
for iRow = 1:size(normDist,1)
    plot(1:3,normDist(iRow,:),'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:3,normDist(iRow,:),'o','color',[0.5 0.5 0.5],'linewidth',1.5)
end
plot(1:3,nanmean(normDist(:,1:3)),'r+','markersize',15,'linewidth',2)
xlim([0.5 3.5])
ylim([0 1])
set(gca,'xtick',[1:3],'xticklabels',{'M1-M2','M1-M3','M2-M3'})
xtickangle(45)
set(gca,'ytick',0:0.2:1)
ylabel('Normalized distance')

%%
load normDistAB_ed5a

figure('position',[723,214,274,314]);
for iRow = 1:size(normDist,1)
    plot(1:3,normDist(iRow,:),'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:3,normDist(iRow,:),'o','color',[0.5 0.5 0.5],'linewidth',1.5)
end
plot(1:3,nanmean(normDist(:,1:3)),'r+','markersize',15,'linewidth',2)
xlim([0.5 3.5])
ylim([0 1])
set(gca,'xtick',[1:3],'xticklabels',{'M1-M2','M1-M3','M2-M3'})
xtickangle(45)
set(gca,'ytick',0:0.2:1)
ylabel('Normalized distance')
