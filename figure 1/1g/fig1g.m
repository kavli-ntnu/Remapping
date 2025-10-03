load expList_MEC_1f
load pkStoreCompiled_1f
load modSpacing_1g

dispStore = cell(1);
normalizedDisp = cell(1);

for iPk = 1:size(pkStoreCompiled,2)
    
    boxSize = expList_MEC{iPk,7}(1);
    numMods = 3;
    pkStore = pkStoreCompiled{1,iPk};
    modSpacing = spacingStore(iPk,1:numMods);
    
    origin = [40 40];
    
    tempDisp = [];
    
    tempDisp(1) = pdist([pkStore(1,1) pkStore(1,2); origin])*boxSize/40;
    tempDisp(2) = pdist([pkStore(2,1) pkStore(2,2); origin])*boxSize/40;
    tempDisp(3) = pdist([pkStore(3,1) pkStore(3,2); origin])*boxSize/40;

    dispStore{iPk,1} = tempDisp;
    
    normalizedDisp{iPk,1} = tempDisp ./ modSpacing(1,1:numMods);
end

%% normalized displacement of each module from origin
normStore = nan(size(normalizedDisp,1),4);

figure('position',[723,214,274,314]);
for iRow = 1:size(normalizedDisp,1)
    thisRow = normalizedDisp{iRow,1};
    
    if size(thisRow,2) == 4
        thisRow = thisRow(1,1:3);
    end
    plot(1:size(thisRow,2),thisRow,'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:size(thisRow,2),thisRow,'o','color',[0.5 0.5 0.5],'linewidth',1.5)
    
    normStore(iRow,1:size(thisRow,2)) = thisRow;
end

xlim([0.5 3.5])
ylim([-0.01 0.6])
set(gca,'xtick',1:3,'xticklabels',{'M1','M2','M3'})
ylabel('Normalized displacement')
box on
plot(1:3,nanmean(normStore(:,1:3)),'r+','markersize',15,'linewidth',2)

%% displacement of each module from origin
dStore = nan(size(dispStore,1),4);

figure('position',[723,214,274,314]);
for iRow = 1:size(dispStore,1)
    thisRow = dispStore{iRow,1};
    
    if size(thisRow,2) == 4
        thisRow = thisRow(1,1:3);
    end
    plot(1:size(thisRow,2),thisRow,'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:size(thisRow,2),thisRow,'o','color',[0.5 0.5 0.5],'linewidth',1.5)
    
    dStore(iRow,1:size(thisRow,2)) = thisRow;
end

xlim([0.5 3.5])
ylim([-2 95])
set(gca,'xtick',1:3,'xticklabels',{'M1','M2','M3'})
set(gca,'ytick',0:10:100)
ylabel('Displacement (cm)')
box on
plot(1:3,nanmean(dStore(:,1:3)),'r+','markersize',15,'linewidth',2)
