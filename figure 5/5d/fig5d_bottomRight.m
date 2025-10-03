load simulationVariables_5d

%%
xVals = [0:0.035:0.5 0.8];
yVals = [0:0.03:0.35 0.5];

var1 = minDist;
var2 = minDisp;
compVar = maxCorrAll;

valMat = [];
for iRow = 1:size(xVals,2)-1
    for iCol = 1:size(yVals,2)-1
        xMin = xVals(iRow);
        xMax = xVals(iRow+1);
        
        yMin = yVals(iCol);
        yMax = yVals(iCol+1);
        
        inds = find(var1 >= xMin & var1 < xMax & var2 >= yMin & var2 < yMax);
        
        if isempty(inds)
            valMat(iRow,iCol) = nan;
        else
            temp = compVar(inds);
            valMat(iRow,iCol) = nanmean(temp);
        end
    end
end

vRot = rot90(valMat);

heatGroup = round(vRot,2);

xvalues = cell(1);
for i = 1:size(xVals,2)-1
    xvalues{1,i} = num2str(xVals(i));
end

yvalues = cell(1);
for i = 1:size(yVals,2)-1
    yvalues{1,i} = num2str(yVals(i));
end
yvalues = fliplr(yvalues);

figure;
heatmap(xvalues,yvalues,heatGroup);

colormap(flipud(jet))
caxis([nanmin(nanmin(heatGroup)) nanmax(nanmax(heatGroup))+0.02])