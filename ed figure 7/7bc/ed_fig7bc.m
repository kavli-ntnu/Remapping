%%
load expList_MEC_1f
load pks1_pks2_ed7b

pkDist = nan(size(pkStore,2),3);
for iExp = 1:size(pkStore,2)
    if ~isempty(pkStore{1,iExp}) && ~isempty(pkStore{2,iExp})
        for iMod = 1:3
            boxSize = expList_MEC{iExp,7}(1);
            nBins = 40;
            pkDist(iExp,iMod) = pdist([pkStore{1,iExp}(iMod,:); pkStore{2,iExp}(iMod,:)])*boxSize/nBins;
        end
    end
end

figure('position',[723,214,274*1.5,314]);
plotSpread({pkDist(:,1) pkDist(:,2) pkDist(:,3)})

xlim([0.5 3.5])
set(gca,'xtick',[1 2 3],'xticklabels',{'M1','M2','M3'})
ylabel('Distance')
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','o','markersize',6,'linewidth',1.5)
hold on
plot(1:3,nanmedian(pkDist(:,1:3)),'r+','markersize',15,'linewidth',2)
box on
ylim([-1 31])

%%
distStoreAB = cell(1);
averageDistAB = nan(size(pkStore,2),1);
for iExp = 1:size(pkStore,2)
    
    boxSize = expList_MEC{iExp,7}(1);
    winners = pkStore{1,iExp}; % for A1B1
    
    tempDist = [];
    
    tempDist(1) = pdist([winners(1,1) winners(1,2); winners(2,1) winners(2,2)])*boxSize/40;
    tempDist(2) = pdist([winners(1,1) winners(1,2); winners(3,1) winners(3,2)])*boxSize/40;
    tempDist(3) = pdist([winners(2,1) winners(2,2); winners(3,1) winners(3,2)])*boxSize/40;

    distStoreAB{iExp,1} = tempDist;
    averageDistAB(iExp,1) = nanmean(tempDist);
end

%%
distStoreA2B2 = cell(1);
averageDistA2B2 = nan(size(pkStore,2),1);
for iExp = 1:size(pkStore,2)
    
    boxSize = expList_MEC{iExp,7}(1);
    winners = pkStore{2,iExp}; % for A2B2
    
    if ~isempty(winners)
        tempDist = [];
        
        tempDist(1) = pdist([winners(1,1) winners(1,2); winners(2,1) winners(2,2)])*boxSize/40;
        tempDist(2) = pdist([winners(1,1) winners(1,2); winners(3,1) winners(3,2)])*boxSize/40;
        tempDist(3) = pdist([winners(2,1) winners(2,2); winners(3,1) winners(3,2)])*boxSize/40;
        
        distStoreA2B2{iExp,1} = tempDist;
        averageDistA2B2(iExp,1) = nanmean(tempDist);
    end
end

%%
figure('position',[-756,272,237,316]);
for iRow = 1:size(averageDistAB,1)
    plot(1:2,[averageDistAB(iRow,1) averageDistA2B2(iRow,1)],'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:2,[averageDistAB(iRow,1) averageDistA2B2(iRow,1)],'o','color',[0.5 0.5 0.5],'linewidth',1.5)
end

xlim([0.5 2.5])
ylim([0 110])
set(gca,'xtick',1:2,'xticklabels',{'A1xB1','A2xB2'})
set(gca,'ytick',0:20:100)
ylabel('Distance (cm)')

%%
load rotation1_rotation2_ed7c

rotDiff = nan(size(rotationA1B1,1),size(rotationA1B1,2));
for iRow = 1:size(rotationA1B1,1)
    for iCol = 1:size(rotationA1B1,2)
        rotDiff(iRow,iCol) = calc.circDiff(rotationA1B1(iRow,iCol),rotationA2B2(iRow,iCol));
    end
end

%% 
figure('position',[723,214,274*1.5,314]);
plotSpread({rotDiff(:,1) rotDiff(:,2) rotDiff(:,3)})

xlim([0.5 3.5])
set(gca,'xtick',[1 2 3],'xticklabels',{'M1','M2','M3'})
ylabel('Rotation difference')
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','o','markersize',6,'linewidth',1.5)
hold on
plot(1:3,nanmedian(rotDiff(:,1:3)),'r+','markersize',15,'linewidth',2)

ylim([-5 105])
set(gca,'ytick',0:20:100)
box on

%%
rotationA2B2(8,3) = rotationA2B2(8,3) + 360;
rotationA1B1(8,3) = rotationA1B1(8,3) + 360;

%%
meanRot1 = nanmean(rotationA1B1,2);
meanRot2 = nanmean(rotationA2B2,2);
expDiff = calc.circDiff(meanRot1,meanRot2);

meanRot2(10,1) = meanRot2(10,1) + 360;

%%
figure('position',[-756,272,237,316]);
for iRow = 1:size(meanRot1,1)
    plot(1:2,[meanRot1(iRow,1) meanRot2(iRow,1)],'-','color',[0.8 0.8 0.8],'linewidth',1.5)
    hold on
    plot(1:2,[meanRot1(iRow,1) meanRot2(iRow,1)],'o','color',[0.5 0.5 0.5],'linewidth',1.5)
end

xlim([0.5 2.5])
ylim([-10 380])
set(gca,'xtick',1:2,'xticklabels',{'A1xB1','A2xB2'})
set(gca,'ytick',0:60:360)
ylabel('Mean rotation (deg)')




