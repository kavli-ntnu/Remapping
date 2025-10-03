load fig3de_sessionPairs
load expList_HP_3de

%%
meanStore = nan(size(expStore,1),1);
errStore = nan(size(expStore,1),1);

for iRow = 1:size(expStore,1)
    if ~isempty(expStore{iRow,5})
        meanStore(iRow,1) = nanmean(expStore{iRow,5});
        errStore(iRow,1) = nanstd(expStore{iRow,5}) ./ sqrt(sum(~isnan(expStore{iRow,5})));
    end
end

meanStoreStab = nan(size(expStore,1),1);
errStoreStab = nan(size(expStore,1),1);

for iRow = 1:size(expStoreStab,1)
    if ~isempty(expStoreStab{iRow,5})
        meanStoreStab(iRow,1) = nanmean(expStoreStab{iRow,5});
        errStoreStab(iRow,1) = nanstd(expStoreStab{iRow,5}) ./ sqrt(sum(~isnan(expStoreStab{iRow,5})));
    end
end

%%
figure;
scatter(minNormDist,meanStore,55,expIdx,'filled')

xlim([-0.02 0.72])
ylim([-0.05 0.8])
set(gca,'ytick',-0.2:0.2:1)
set(gca,'xtick',0:0.1:0.7)
xlabel('Minimum distance')
ylabel('Spatial correlation')
box on

