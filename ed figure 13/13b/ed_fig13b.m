load displacement_remapping

bins = 0:0.05:0.6;

store = nan(1,size(bins,2));
numStore = nan(1,size(bins,2));
stdStore = nan(1,size(bins,2));

testVar = minDisp;
d = testVar;
compVar = maxCorrAll;


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

xVals = bins(1:size(bins,2));
yVals = store(1:size(bins,2));

createFit13BL(xVals,yVals)
set(gcf,'color','w')
set(gcf,'position',[-710,144,560,420])
grid off
legend off
xlim([-0.02 0.52])
ylim([0.2 0.8])
set(gca,'xtick',0:0.1:0.5,'ytick',0.2:0.1:1)
xlabel('Minimum displacement')
ylabel('PV correlation')


for i = 1:size(bins,2)
    dMin = bins(i);
    
    if i == size(bins,2)
        temp = maxCorrAll(d >= dMin);
    else
        dMax = bins(i+1);
        temp = maxCorrAll(d >= dMin & d < dMax);
    end

    store(i) = nanmean(temp);

    hold on
    errorbar(xVals(i),nanmean(store(i)),(nanstd(temp)) ./ (sqrt(sum(~isnan(temp)))),'-','markersize',10,'markeredgecolor','k','markerfacecolor','k','linewidth',1.5,'color','k');
    hold on
end

%%
bins = 0:0.05:0.6;

store = nan(1,size(bins,2));
numStore = nan(1,size(bins,2));
stdStore = nan(1,size(bins,2));

testVar = minDisp;
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

xVals = bins(1:size(bins,2));
yVals = store(1:size(bins,2));

createFit13BM(xVals,yVals)
set(gcf,'color','w')
set(gcf,'position',[-710,144,560,420])
grid off
legend off
xlim([-0.02 0.52])
ylim([0.6 1])
set(gca,'xtick',0:0.1:0.5,'ytick',0.2:0.1:1)
xlabel('Minimum displacement')
ylabel('Rearrangment score')


for i = 1:size(bins,2)
    dMin = bins(i);
    
    if i == size(bins,2)
        temp = rValAll(d >= dMin);
    else
        dMax = bins(i+1);
        temp = rValAll(d >= dMin & d < dMax);
    end

    store(i) = nanmean(temp);

    hold on
    errorbar(xVals(i),nanmean(store(i)),(nanstd(temp)) ./ (sqrt(sum(~isnan(temp)))),'-','markersize',10,'markeredgecolor','k','markerfacecolor','k','linewidth',1.5,'color','k');
    hold on
end

%%
bins = 0:0.05:0.6;

store = nan(1,size(bins,2));
numStore = nan(1,size(bins,2));
stdStore = nan(1,size(bins,2));

testVar = minDisp;
d = testVar;
compVar = validStoreAll;

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

xVals = bins(1:size(bins,2));
yVals = store(1:size(bins,2));

%%
createFit13BR(xVals,yVals)
set(gcf,'color','w')
set(gcf,'position',[-710,144,560,420])
grid off
legend off
xlim([-0.02 0.52])
ylim([0.32 0.59])
set(gca,'xtick',0:0.1:0.5,'ytick',0.2:0.1:1)
xlabel('Minimum displacement')
ylabel('Percent turnover')

for i = 1:size(bins,2)
    dMin = bins(i);
    
    if i == size(bins,2)
        temp = validStoreAll(d >= dMin);
    else
        dMax = bins(i+1);
        temp = validStoreAll(d >= dMin & d < dMax);
    end

    store(i) = nanmean(temp);

    hold on
    errorbar(xVals(i),nanmean(store(i)),(nanstd(temp)) ./ (sqrt(sum(~isnan(temp)))),'-','markersize',10,'markeredgecolor','k','markerfacecolor','k','linewidth',1.5,'color','k');
    hold on
end
