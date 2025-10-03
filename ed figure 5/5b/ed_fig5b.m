load distStore_1h
dAll = distStore;

load spacingAllExpt_ed5a

theoreticalMax = nan(size(spacingAllExpt,1),3);
for iExp = 1:size(spacingAllExpt,1)
    theoreticalMax(iExp,1) = spacingAllExpt(iExp,2) * tan(pi/6);
    theoreticalMax(iExp,2:3) = spacingAllExpt(iExp,3) * tan(pi/6);
end

load finalDistStore_1h

expPrc = nan(size(finalDistStore,1),3);
for iRow = 1:size(finalDistStore,1)
    if ~isempty(finalDistStore{iRow,1})
        f = finalDistStore{iRow,1};
        
        expPrc(iRow,:) = prctile(f,5);
    end
end
eAll = expPrc;

expMedian = nan(size(finalDistStore,1),3);
for iRow = 1:size(finalDistStore,1)
    if ~isempty(finalDistStore{iRow,1})
        f = finalDistStore{iRow,1};
        
        expMedian(iRow,:) = nanmedian(f);
    end
end

normPrc = nan(size(finalDistStore,1),3);
normMedian = nan(size(finalDistStore,1),3);
for iRow = 1:size(finalDistStore,1)
    if ~isempty(finalDistStore{iRow,1})
        temp = finalDistStore{iRow,1};    
        tempNorm = temp ./ theoreticalMax(iRow,:);
        
        normPrc(iRow,:) = prctile(tempNorm,5);
        normMedian(iRow,:) = nanmedian(tempNorm);
    end
end

%%
diff5 = dAll - eAll;

figure('position',[659,157,457,357]);

plotSpread(diff5)
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','o','markersize',6,'linewidth',1.5)
xlim([0.5 3.5])
ylim([-3 70])
yline(0)
box on
set(gca,'ytick',0:15:100,'xtick',1:3,'xticklabels',{'M1-M2','M1-M3','M2-M3'})
ylabel('Difference')

%%
diff50 = dAll - expMedian;

figure('position',[659,157,457,357]);
plotSpread(diff50)
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','o','markersize',6,'linewidth',1.5)
xlim([0.5 3.5])
ylim([-40 40])
yline(0)
box on
set(gca,'ytick',-40:20:100,'xtick',1:3,'xticklabels',{'M1-M2','M1-M3','M2-M3'})
ylabel('Difference')

%% 
load normDistAB_ed5a

diff5 = normDist - normPrc;

figure('position',[-659,157,457,357]);
plotSpread(diff5)
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','o','markersize',6,'linewidth',1.5)
xlim([0.5 3.5])
ylim([-0.05 0.8])
yline(0)
box on
set(gca,'ytick',0:0.2:1,'xtick',1:3,'xticklabels',{'M1-M2','M1-M3','M2-M3'})
ylabel('Difference')

%%
diff50 = normDist - normMedian;

figure('position',[659,157,457,357]);
plotSpread(diff50)
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','o','markersize',6,'linewidth',1.5)
xlim([0.5 3.5])
ylim([-0.45 0.45])
yline(0)
box on
set(gca,'ytick',-0.4:0.2:1,'xtick',1:3,'xticklabels',{'M1-M2','M1-M3','M2-M3'})
ylabel('Difference')

