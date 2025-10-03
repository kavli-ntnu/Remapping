load fig3de_sessionPairs
load expList_HP_3de

%%
figure('position',[723,214,274,314]);

for iExp = 1:size(normDist,1)  
    temp = normDist(iExp,:);
    
    plot(1:size(temp,2),temp,'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:size(temp,2),temp,'o','color',[0.5 0.5 0.5],'linewidth',1.5)   
end

xlim([0.5 3.5])
ylim([-0.02 1.02])
set(gca,'xtick',1:3,'xticklabels',{'M1-M2','M1-M3','M2-M3'})
set(gca,'ytick',0:0.2:2)
ylabel('Normalized distance')
box on
