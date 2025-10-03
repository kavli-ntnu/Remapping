load rotation_4cd

%%
rotFig = rotation;
rotFig(5,3) = rotFig(5,3) + 360;
rotFig(10,3) = rotFig(10,3) + 360;
rotFig(14,1:2) = rotFig(14,1:2) + 360;
rotFig(15,2) = rotFig(15,2) + 360;
rotFig(19,3) = rotFig(19,3) + 360;
rotFig(27,3) = rotFig(27,3) + 360;

%%
figure('position',[-1858,497,274,314]);
for iRow = 1:size(rotFig,1)
    plot(1:3,rotFig(iRow,1:3),'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:3,rotFig(iRow,1:3),'o','color',[0.5 0.5 0.5],'linewidth',1.5)
end

xlim([0.5 3.5])
ylim([-10 380])
set(gca,'ytick',0:60:360)
set(gca,'xtick',1:3,'xticklabels',{'M1','M2','M3'})
ylabel('Module rotation (deg)')

%%
column1_all = calc.circDiff(rotation(:,2),rotation(:,1));
column2_all = calc.circDiff(rotation(:,3),rotation(:,1));
column3_all = calc.circDiff(rotation(:,3),rotation(:,2));
 
allDiff = [column1_all column2_all column3_all];
maxDiff = nanmax(allDiff,[],2);

allCol = vertcat(column1_all,column2_all,column3_all);

nanmedian(allCol)
sum(~isnan(allCol))

%%
figure('position',[-1858,80,411,314]);
plotSpread({column1_all column2_all column3_all})

xlim([0.5 3.5])
set(gca,'xtick',[1 2 3],'xticklabels',{'M1-M2','M1-M3','M2-M3'})
ylabel('Rotation difference')
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','.','markersize',10,'linewidth',1.5)
hold on
plot(1:3,[nanmedian(column1_all) nanmedian(column2_all) nanmedian(column3_all)],'r+','markersize',15,'linewidth',2)

nanmedian(allCol)
sum(~isnan(allCol))

ylim([-2 92])
set(gca,'ytick',0:20:90)
plot([0 90],[90 90],'--','color',[0.8 0.8 0.8])
box on

