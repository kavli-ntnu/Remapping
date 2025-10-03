load simulationVariables_5d

%%
bins = [0:0.038:0.152 0.21 0.266:0.038:1]; % minDist

store = nan(1,size(bins,2));
numStore = nan(1,size(bins,2));
stdStore = nan(1,size(bins,2));

testVar = minDist;
d = testVar;
compVar = rValAll;

for i = 1:size(bins,2)
    dMin = bins(i);
    
    if i == size(bins,2)
        temp = compVar(d >= dMin);
    else
        dMax = bins(i+1);
        temp = compVar(d >= dMin & d < dMax);
    end
    
    store(i) = nanmean(temp);
    numStore(i) = size(temp,1);
    stdStore(i) = nanstd(temp);
end

%%
xVals = bins(1:size(bins,2));
yVals = store(1:size(bins,2));
createFit1(xVals,yVals)

set(gcf,'color','w')
grid off
legend off
xlim([-0.02 0.8])
ylim([0.6 1])
set(gca,'xtick',0:0.2:1,'ytick',0.2:0.1:1)
xlabel('Minimum distance')
ylabel('Rearrangement score')

%%
for i = 1:size(bins,2)
    dMin = bins(i);
    
    if i == size(bins,2)
        temp = compVar(d >= dMin);
    else
        dMax = bins(i+1);
        temp = compVar(d >= dMin & d < dMax);
    end

    store(i) = nanmean(temp);

    hold on
    errorbar(xVals(i),nanmean(store(i)),(nanstd(temp)) ./ (sqrt(sum(~isnan(temp)))),'-','markersize',10,'markeredgecolor','k','markerfacecolor','k','linewidth',1.5,'color','k');
    hold on
end

