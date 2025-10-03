load m4_ed11j

figure('position',[858,497,274,314]);
for iRow = 1:4
    plot(1:4,m4(iRow,1:4),'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:4,m4(iRow,1:4),'o','color',[0.5 0.5 0.5],'linewidth',1.5)
end

xlim([0.5 4.5])
ylim([-10 380])
set(gca,'ytick',0:60:360)
set(gca,'xtick',1:4,'xticklabels',{'M1','M2','M3','M4'})
ylabel('Module rotation (deg)')
