load rotationAA_ed11a

%%
column1 = calc.circDiff(rotation(:,2),rotation(:,1));
column2 = calc.circDiff(rotation(:,3),rotation(:,1));
column3 = calc.circDiff(rotation(:,3),rotation(:,2));

allDiff = [column1 column2 column3];

figure; 
boxplot(allDiff)
ylim([-2 92])
set(gca,'ytick',0:20:90)
hold on
plot([0 90],[90 90],'--','color',[0.8 0.8 0.8])
box on
ylabel('Rotation difference')
set(gca,'xtick',1:3,'xticklabels',{'M1-M2','M1-M3','M2-M3'})

%%
rotation(10,3) = rotation(10,3) - 360; % rotations extend across 360

%%
figure('position',[723,214,274,314]);
for iRow = 1:size(rotation,1)
    plot(1:3,rotation(iRow,1:3),'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:3,rotation(iRow,1:3),'o','color',[0.5 0.5 0.5],'linewidth',1.5)
end

xlim([0.5 3.5])
ylim([-10 380])
set(gca,'xtick',1:3,'xticklabels',{'M1','M2','M3'})
set(gca,'ytick',0:60:360)
ylabel('Module rotation (deg)')
