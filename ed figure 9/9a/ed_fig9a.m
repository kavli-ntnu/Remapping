%%
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

%%
gridVar = minNormDist;
secondVar = meanStore;
X = [gridVar secondVar];

maxK = 7;
meanS = nan(maxK,1);

for k = 1:maxK
    [idx,C] = kmeans(X,k,'distance','cityblock','replicates',1000);

    load kMeansColorVec
    colorMat = colorVec/255;
    
    subplot(121)
    for c = 1:k
        inds = find(idx == c);
        scatter(gridVar(inds),secondVar(inds),55,'filled')
        hold on
        scatter(C(c,1),C(c,2),'k*')
    end
    
    ptDist = cell(1);
    for c = 1:k
        inds = find(idx == c);
        
        tempDist = nan(size(inds,1),1);
        for i = 1:size(inds,1)
            tempDist(i,1) = pdist([C(c,1) C(c,2); gridVar(inds(i)) secondVar(inds(i))],'cityblock');
        end
        
        sqTemp = tempDist .* tempDist;
        ptDist{c,1} = sqTemp;
    end
    ptMat = cell2mat(ptDist);
    sse(k,1) = nansum(ptMat);
        
    % X = X(~isnan(X(:,1)),:);
    % idx = idx(~isnan(idx));
    
    subplot(122)
    [s,h] = silhouette(X,idx,'cityblock');
    close(gcf)
    allS = nanmean(s);
    allSvals(k,1) = allS;
    
    clusterS = [];
    for c = 1:k
        clusterS(c,1) = nanmean(s(idx == c));
    end
    
    % clusterS
    meanS(k,1) = nanmean(clusterS);

    %%
    figure;
    for i = 1:k
        currentIdx = find(idx == i);
        
        corrCol = [];
        for iRow = 1:size(currentIdx,1)
            corrCol = [corrCol; expStore{currentIdx(iRow),5}];
        end
        
        [f,xi] = ksdensity(corrCol,'bandwidth',0.11);
        [pks,loc,width,prom] = findpeaks(f,xi);
        findpeaks(f,xi)
        hold on
    end
    close(gcf)
end

figure;
plot(1:size(sse,1),sse,'-o')
xlim([0.5 maxK+0.5])
xlabel('k')
ylabel('sse')
set(gca,'ytick',0:1:10)