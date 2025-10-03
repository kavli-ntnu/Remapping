load table3_grpStoreMEC

medianStore = nan(size(grpStore,1),size(grpStore,2));
for iRow = 1:size(grpStore,1)
    for iGrp = 1:size(grpStore,2)
        temp = grpStore{iRow,iGrp};
        medianStore(iRow,iGrp) = nanmedian(temp);
    end
end

%%
nPerGroup = nan(size(grpStore,1),size(grpStore,2));
for iRow = 1:size(grpStore,1)
    for iGrp = 1:size(grpStore,2)
        temp = grpStore{iRow,iGrp};
        nPerGroup(iRow,iGrp) = sum(~isnan(temp));
    end
end
nanmax(nPerGroup)