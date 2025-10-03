load rotShift_ed13e

%%
compVar = maxCorr;

vec = [0.8:-0.1:0]';
colorMap = [vec vec vec];
rotGrp = 1:1500:13500;

figure;
bins = [0.01:0.06:0.5 1];

store = nan(1,size(bins,2));
numStore = nan(1,size(bins,2));
stdStore = nan(1,size(bins,2));
finalStore = [];

for iGrp = 1:size(rotGrp,2)
    inds = rotGrp(iGrp):rotGrp(iGrp)+1499;
    
    testVar = minDist(inds,1);
    d = testVar;
    
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
    
    finalStore(iGrp,:) = store;
    
    xVals = bins(1:size(bins,2));
    yVals = store(1:size(bins,2));
    
    xLogic = ~isnan(xVals);
    yLogic = ~isnan(yVals);
    xVals = xVals(xLogic & yLogic);
    yVals = yVals(xLogic & yLogic);

    f = fit(xVals(:), yVals(:), 'smoothingspline');
    h = plot(f);
    set(h(1),'color',colorMap(iGrp,:),'linewidth',1.5); 
    
    hold on
    plot(xVals,yVals,'.','color',colorMap(iGrp,:),'markersize',25)
    legend off
end

ylim([0.2 0.8])
xlabel('Minimum distance')
ylabel('PV correlation')
set(gca,'xtick',0:0.1:1,'ytick',0.2:0.1:1)
xlim([-0.02 0.52])

%%
compVar = remappingVal;

vec = [0.8:-0.1:0]';
colorMap = [vec vec vec];
rotGrp = 1:1500:13500;

figure;
bins = [0.01:0.06:0.5 1];

store = nan(1,size(bins,2));
numStore = nan(1,size(bins,2));
stdStore = nan(1,size(bins,2));
finalStore = [];

for iGrp = 1:size(rotGrp,2)
    inds = rotGrp(iGrp):rotGrp(iGrp)+1499;
    
    testVar = minDist(inds,1);
    d = testVar;
    
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
    
    finalStore(iGrp,:) = store;
    
    xVals = bins(1:size(bins,2));
    yVals = store(1:size(bins,2));
    
    xLogic = ~isnan(xVals);
    yLogic = ~isnan(yVals);
    xVals = xVals(xLogic & yLogic);
    yVals = yVals(xLogic & yLogic);

    f = fit(xVals(:), yVals(:), 'smoothingspline');
    h = plot(f);
    set(h(1),'color',colorMap(iGrp,:),'linewidth',1.5); 
    
    hold on
    plot(xVals,yVals,'.','color',colorMap(iGrp,:),'markersize',25)
    legend off
end

ylim([0.6 1])
xlabel('Minimum distance')
ylabel('Remapping strength')
set(gca,'xtick',0:0.1:1,'ytick',0.6:0.1:1)
xlim([-0.02 0.52])

%%
compVar = validStoreAll;

vec = [0.8:-0.1:0]';
colorMap = [vec vec vec];
rotGrp = 1:1500:13500;

figure;
bins = [0.01:0.06:0.5 1];

store = nan(1,size(bins,2));
numStore = nan(1,size(bins,2));
stdStore = nan(1,size(bins,2));
finalStore = [];

for iGrp = 1:size(rotGrp,2)
    inds = rotGrp(iGrp):rotGrp(iGrp)+1499;
    
    testVar = minDist(inds,1);
    d = testVar;
    
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
    
    finalStore(iGrp,:) = store;
    
    xVals = bins(1:size(bins,2));
    yVals = store(1:size(bins,2));
    
    xLogic = ~isnan(xVals);
    yLogic = ~isnan(yVals);
    xVals = xVals(xLogic & yLogic);
    yVals = yVals(xLogic & yLogic);

    f = fit(xVals(:), yVals(:), 'smoothingspline');
    h = plot(f);
    set(h(1),'color',colorMap(iGrp,:),'linewidth',1.5); 
    
    hold on
    plot(xVals,yVals,'.','color',colorMap(iGrp,:),'markersize',25)
    legend off
end

ylim([0.35 0.55])
xlabel('Minimum distance')
ylabel('Percent turnover')
set(gca,'xtick',0:0.1:1,'ytick',0.35:0.05:1)
xlim([-0.02 0.52])
