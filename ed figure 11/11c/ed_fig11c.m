load moduleOrient_ed11c

%%
figure;
plotSpread(abs(meanDiffA))
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','o','markersize',6,'linewidth',1.5)
ylim([-1 30])
hold on
box on
plot(1:3,nanmedian(abs(meanDiffA)),'r+','markersize',15,'linewidth',1.5)
set(gca,'xticklabels',{'M1-M2','M1-M3','M2-M3'})
ylabel('Orientation difference')

%%
diffAB(:,1) = meanDiffA(:,1)-meanDiffB(:,1);
diffAB(:,2) = meanDiffA(:,2)-meanDiffB(:,2);
diffAB(:,3) = meanDiffA(:,3)-meanDiffB(:,3);

figure;
plotSpread(abs(diffAB))
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','o','markersize',6,'linewidth',1.5)
ylim([-1 30])
hold on
box on
plot(1:3,nanmedian(abs(diffAB)),'r+','markersize',15,'linewidth',1.5)
set(gca,'xticklabels',{'M1-M2','M1-M3','M2-M3'})
xlim([0.25 3.75])
ylabel('Change in offset')
title('A1xB1')

%%
diffAA(:,1) = meanDiffA(:,1)-meanDiffA2(:,1);
diffAA(:,2) = meanDiffA(:,2)-meanDiffA2(:,2);
diffAA(:,3) = meanDiffA(:,3)-meanDiffA2(:,3);

figure;
plotSpread(abs(diffAA))
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','o','markersize',6,'linewidth',1.5)
ylim([-1 30])
hold on
box on
plot(1:3,nanmedian(abs(diffAA)),'r+','markersize',15,'linewidth',1.5)
set(gca,'xticklabels',{'M1-M2','M1-M3','M2-M3'})
xlim([0.25 3.75])
ylabel('Change in offset')
title('A1xA2')
