
load table2_fxStore
% row 1 = grid
% row 2 = border
% row 3 = HD
% row 4 = spatial
% row 5 = other

% columns listed in header variable

%% mean of one measure for one group
currentCol = 8

for iRow = 1 %:size(fxStore,1)
    currentMat = [];
    if ~isempty(fxStore{iRow,currentCol})
        temp = fxStore{iRow,currentCol};

        if size(temp,2) == 4
            currentMat = [currentMat; temp];
        elseif size(temp,2) == 2
            tempMat = nan(size(temp,1),4);
            tempMat(:,1) = temp(:,1);
            tempMat(:,3) = temp(:,2);
            currentMat = [currentMat; tempMat];
        elseif size(temp,2) == 1
            tempMat = nan(size(temp,1),4);     
            tempMat(:,3) = temp;
            currentMat = [currentMat; tempMat]; 
        end 
    end
end

if currentCol == 3 || currentCol == 4 || currentCol == 9
    currentMat = abs(currentMat);
end

nanmedian(currentMat)

%%
for iCol = 1:size(currentMat,2)
    data = currentMat(:,iCol);
    data = data(~isnan(data));
    nBoot = 1000;
    
    bootMedians = zeros(nBoot,1);
    for i = 1:nBoot
        resampleData = datasample(data,length(data));
        bootMedians(i) = median(resampleData);
    end
    ci = prctile(bootMedians,[2.5 97.5])
    med = nanmedian(data);    
end

%%
[p,~,stats] = ranksum(currentMat(:,1),currentMat(:,3))
d = computeCohen_d(currentMat(:,1),currentMat(:,3))

