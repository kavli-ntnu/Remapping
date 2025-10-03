load prcTurnover_ed8b

figure('position',[1285,67,198,420]);
plotSpread({prcStab prcRemap})
set(findobj(gca,'type','line'),'color','k','marker','o','markersize',6,'linewidth',1.5)
hold on
plot(1:2,[nanmean(prcStab) nanmean(prcRemap)],'r+','markersize',15,'linewidth',2)
set(gca,'xtick',1:2,'ytick',0.1:0.1:0.5)
xlim([0.5 2.5])
ylim([0.1 0.5])
ylabel('Percent turnover')
