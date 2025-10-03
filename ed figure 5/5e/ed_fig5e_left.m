load moduleSeparation_ed5e

%%
sepStore = nan(size(depthLocAll,1),3);
for iRow = 1:size(depthLocAll,1)
    if ~isempty(depthLocAll{iRow,1})
        d = depthLocAll{iRow,1};
        s = spacingLocAll{iRow,1};
        
        coords = [d s];
        
        sepStore(iRow,1) = pdist([coords(1,:); coords(2,:)]);
        sepStore(iRow,2) = pdist([coords(1,:); coords(3,:)]);
        sepStore(iRow,3) = pdist([coords(2,:); coords(3,:)]);
    end
end

%%
load distStore_1h

dCol = reshape(distStore,[],1);
sCol = reshape(sepStore,[],1);

figure;
scatter(sCol,dCol)
xlabel('Module separation')
ylabel('Distance between modules')
box on
xlim([0 1250])
ylim([0 90])
set(gca,'ytick',0:20:100)

xVals = sCol;
yVals = dCol;
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);
[r_val p_val] = corr(xVals,yVals)

