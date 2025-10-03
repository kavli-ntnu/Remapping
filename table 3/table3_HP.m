load table3_grpStoreHP

%%
medianStore = nan(size(pcStore,1),size(pcStore,2));
for iRow = 1:size(pcStore,1)
    for iGrp = 1:size(pcStore,2)
        if size(pcStore{iRow,iGrp},2) == 1
            temp = pcStore{iRow,iGrp};
        else
            temp = pcStore{iRow,iGrp}(:,2);
        end
        medianStore(iRow,iGrp) = nanmedian(temp);
    end
end

%%
nPerGroup = nan(size(pcStore,1),size(pcStore,2));
for iRow = 1:size(pcStore,1)
    for iGrp = 1:size(pcStore,2)
        temp = pcStore{iRow,iGrp}(:,1);
        nPerGroup(iRow,iGrp) = sum(~isnan(temp));
    end
end
nanmax(nPerGroup)
