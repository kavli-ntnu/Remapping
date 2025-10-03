load rearrScore_ed8b

figure('position',[1285,67,138,420]);
plotSpread({stab remap})
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','o','markersize',6,'linewidth',1.5)
hold on
plot(1:2,[nanmean(stab) nanmean(remap)],'r+','markersize',15,'linewidth',2)
set(gca,'ytick',0.4:0.2:1.2)
% ylim([0.35 1.1])
ylabel('Rearrangement score')