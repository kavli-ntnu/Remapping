load acStore_26821
boxSize = 200;

numMods = 3;
figure;
spacingStore = nan(1,numMods);

for iMod = 1:numMods
    subplot(1,numMods,iMod)
    colorMapBRK(rotCCstore(:,:,iMod));
    [~,gStats] = analyses.gridnessScore(rotCCstore(:,:,iMod),'threshold',0.1);
    szmap = (size(rotCCstore,1)/2);
    spacingStore(1,iMod) = nanmean(gStats.spacing)*boxSize/szmap;
    hold on
    COM = comStore{1,iMod};
    plot(COM(:,1),COM(:,2),'wo')
end
spacingStore

%%
load acStore_27150
boxSize = 150;

numMods = 4;
spacingStore = nan(1,numMods);

figure;
for iMod = 1:numMods
    subplot(1,numMods,iMod)
    colorMapBRK(rotCCstore(:,:,iMod));
    [~,gStats] = analyses.gridnessScore(rotCCstore(:,:,iMod),'threshold',0.1);
    szmap = (size(rotCCstore,1)/2);
    spacingStore(1,iMod) = nanmean(gStats.spacing)*boxSize/szmap;
    hold on
    COM = comStore{1,iMod};
    if ~isempty(COM)
        plot(COM(:,1),COM(:,2),'wo')
    end
end
spacingStore
