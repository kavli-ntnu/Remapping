load fig3de_sessionPairs
load expList_HP_3de
grpIDs = [1 3 2];

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

meanSC = meanStore;
meanStabSC = meanStoreStab;

expStoreSC = expStore;
expStoreStabSC = expStoreStab;

%%
load expListMEC_ed9
load expStoreMEC_ed9

colString = {'Exp','Num','Mean RD','Peak RD','SI diff','Size diff','Change in grid spacing','Ellip diff','Change in grid score'};
 
figure('position',[-1280,-34,1196,923]);
selectedCols = [3:9];

statStore = nan(size(selectedCols,2),2);

for i = 1:size(selectedCols,2)
    col = selectedCols(i);
    
    colStore = cell(1);
    for iExp = 1:size(expStore,1)
        if ~isempty(expStore{iExp,col})
            expCol = expList_MEC{iExp,8};
            colStore{iExp,1} = expStore{iExp,col}(:,expCol);
        end
    end
    
    rateMatFinal = [];
    grpSCfinal = [];
    
    for iGrp = grpIDs
        rateMat = nan(size(colStore,1),1);
        inds = find(idx == grpIDs(iGrp));
        
        for j = 1:size(inds,1)
            if col == 3 || col == 4
                temp = abs(colStore{inds(j),1});
            else
                temp = colStore{inds(j),1};
            end
            rateMat(j,1) = nanmean(temp);
        end
        
        rateMat = rateMat(1:size(inds,1),1);
        grpSC = meanSC(inds,1);
        
        subplot(3,3,i)
        scatter(rateMat,grpSC,35)
        xlabel(sprintf('%s',colString{1,col}))
        ylabel('Spatial correlation')
        box on
        hold on
        
        rateMatFinal = [rateMatFinal; rateMat];
        grpSCfinal = [grpSCfinal; grpSC];
    end

    xVals = rateMatFinal;
    yVals = grpSCfinal;
    xLogic = ~isnan(xVals);
    yLogic = ~isnan(yVals);
    xVals = xVals(xLogic & yLogic);
    yVals = yVals(xLogic & yLogic);
    
    [r,p] = corr(xVals,yVals);
    statStore(i,1) = r;
    statStore(i,2) = p;

    if col == 3
        subplot(331)
        xlim([0 0.4])
        ylim([-0.02 0.8])
    elseif col == 4
        subplot(332)
        xlim([0 0.4])
        ylim([-0.02 0.8])
    elseif col == 5
        subplot(333)
        xlim([-0.2 0.2])
        ylim([-0.02 0.8])
    elseif col == 6
        subplot(334)
        xlim([-0.2 0.2])
        ylim([-0.02 0.8])
    elseif col == 7
        subplot(335)
        xlim([-0.15 0.15])
        ylim([-0.02 0.8])
    elseif col == 8
        subplot(336)
        xlim([-0.62 0.62])
        ylim([-0.02 0.8])
    elseif col == 9
        subplot(337)
        xlim([-0.52 0.52])
        ylim([-0.02 0.8])
    end
end

