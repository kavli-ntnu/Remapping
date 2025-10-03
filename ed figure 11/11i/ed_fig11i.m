load modRotation_ed11i
groupIDs = [2 3 1];

figure('position',[1490,302,1476,358]);
subplot(141)
rot1 = rotFig(idx == groupIDs(1),:);
for iRow = 1:size(rot1,1)
    plot(1:3,rot1(iRow,1:3),'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:3,rot1(iRow,1:3),'o','color',[0.5 0.5 0.5],'linewidth',1.5)
end
xlim([0.5 3.5])
ylim([-10 380])
set(gca,'ytick',0:60:360)
set(gca,'xtick',1:3,'xticklabels',{'M1','M2','M3'})
ylabel('Module rotation (deg)')

subplot(142)
rot2 = rotFig(idx == groupIDs(2),:);
for iRow = 1:size(rot2,1)
    plot(1:3,rot2(iRow,1:3),'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:3,rot2(iRow,1:3),'o','color',[0.5 0.5 0.5],'linewidth',1.5)
end
xlim([0.5 3.5])
ylim([-10 380])
set(gca,'ytick',0:60:360)
set(gca,'xtick',1:3,'xticklabels',{'M1','M2','M3'})

subplot(143)
rot3 = rotFig(idx == groupIDs(3),:);
for iRow = 1:size(rot3,1)
    plot(1:3,rot3(iRow,1:3),'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:3,rot3(iRow,1:3),'o','color',[0.5 0.5 0.5],'linewidth',1.5)
end
xlim([0.5 3.5])
ylim([-10 380])
set(gca,'ytick',0:60:360)
set(gca,'xtick',1:3,'xticklabels',{'M1','M2','M3'})

subplot(144)
for iRow = 1:size(rotFigAA,1)
    plot(1:3,rotFigAA(iRow,1:3),'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:3,rotFigAA(iRow,1:3),'o','color',[0.5 0.5 0.5],'linewidth',1.5)
end

xlim([0.5 3.5])
ylim([-10 380])
set(gca,'ytick',0:60:360)
set(gca,'xtick',1:3,'xticklabels',{'M1','M2','M3'})
ylabel('Module rotation (deg)')

%%
figure;
plotSpread({grpRotation{1,1} grpRotation{1,2} grpRotation{1,3}})
hold on
plot(1:3,[nanmedian(grpRotation{1,1}) nanmedian(grpRotation{1,2}) nanmedian(grpRotation{1,3})],'r+','markersize',15,'linewidth',2)
ylim([-2 92])
box on
ylabel('Rotation difference')
