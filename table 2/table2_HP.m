load table2_pcStore

%% columns: 1 = all, 2 = CA3, 3 = CA1, 4 = DG
selectedCol = 1;

rowMedians = nan(size(pcStore,1),2);
for iRow = 1:size(pcStore,1)
    temp = pcStore{iRow,selectedCol};
    
    rowMedians(iRow,:) = [nanmedian(temp{1,1}) nanmedian(temp{2,1})]; 
end
    