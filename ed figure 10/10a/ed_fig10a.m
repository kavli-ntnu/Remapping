load fig3de_sessionPairs
load expList_HP_3de
groupIDs = [1 3 2];

load colorMat

%%
xStore = cell(1);
yStore = cell(1);
countStore = cell(1);
nBins = 100;

figure;
for iGrp = 1:3
    inds = find(idx == groupIDs(iGrp));
    
    x1 = nan(size(inds,1),nBins);
    y1 = nan(size(inds,1),nBins);
    counts1 = nan(size(inds,1),nBins);
    
    for i = 1:size(inds,1)
        temp = expStore{inds(i),5};
        h = histfit(temp,nBins,'kernel');
        x1(i,:) = h(2).XData;
        y1(i,:) = h(2).YData;
        counts1(i,:) = histcounts(temp,nBins);
    end
    
    xStore{1,iGrp} = x1;
    yStore{1,iGrp} = y1;
    countStore{1,iGrp} = counts1;
end
close(gcf)

figure('position',[-1719,38,1663,297]);
for iGrp = 1:3
    subplot(1,4,iGrp)
    
    xTemp = xStore{1,iGrp};
    yTemp = yStore{1,iGrp};
    countTemp = countStore{1,iGrp};
    
    for i = 1:size(xTemp,1)
        x = xTemp(i,:);
        y = yTemp(i,:);
        c = countTemp(i,:);
        
        plot(xTemp(i,:),y/sum(c),'-','color',colorMat(iGrp,:))
        hold on
    end
    ylim([0 0.03])
    xlim([-0.5 1.05])
    ylabel('Frequency')
    xlabel('Spatial correlation')
end

% stable
stabCorrStore = cell(1);
stabCorr = [];
for iRow = 1:size(expStoreStab,1)
    if ~isempty(expStoreStab{iRow,5}) % stab col
        stabCorrStore{iRow,1} = expStoreStab{iRow,5};
        stabCorr = [stabCorr; expStoreStab{iRow,5}];
    end
end

x1 = nan(size(stabCorrStore,1),nBins);
y1 = nan(size(stabCorrStore,1),nBins);
counts1 = nan(size(stabCorrStore,1),nBins);

for iRow = 1:size(stabCorrStore,1)
    if ~isempty(stabCorrStore{iRow,1})
        temp = stabCorrStore{iRow,1};
        
        if sum(~isnan(temp))
            figure;
            h = histfit(temp,nBins,'kernel');
            x1(iRow,:) = h(2).XData;
            y1(iRow,:) = h(2).YData;
            counts1(iRow,:) = histcounts(temp,nBins);
            close(gcf)
        end
    end
end

x1 = minions.removeNans(x1,'rows','all');
y1 = minions.removeNans(y1,'rows','all');
counts1 = minions.removeNans(counts1,'rows','all');

subplot(144)
for i = 1:size(x1,1)
    x = x1(i,:);
    y = y1(i,:);
    c = counts1(i,:);
    
    plot(x,y/sum(c),'-','color',[0.5 0.5 0.5])
    hold on
    
    stabX = x;
    stabY = y;
    stabC = c;
end
ylim([0 0.03])
xlim([-0.5 1.05])
ylabel('Frequency')
xlabel('Spatial correlation')