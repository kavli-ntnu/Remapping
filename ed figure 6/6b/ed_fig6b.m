load rotDiff_ed6b

figure('position',[-1778,528,1311,262]);
for iMod = 1:3
    rotDiffAll = nan(size(rotDiffFinal,1),25);
    for iRow = 1:size(rotDiffFinal,1)
        temp = rotDiffFinal{iRow,iMod};
        rotDiffAll(iRow,1:size(temp,1)) = temp;
    end
    rotDiffAll = abs(rotDiffAll);
     
    cDiff = nan(size(rotDiffAll,1),size(rotDiffAll,2));
    for iCol = 1:size(rotDiffAll,2)
        temp = rotDiffAll(:,iCol);
        testCol = repmat(360,[size(temp,1) 1]);
        cDiff(:,iCol) = calc.circDiff(temp,testCol,1);
    end
    
    xVals = 5:5:size(cDiff,2)*5;
    for iRow = 1:size(cDiff,1)
        yVals = cDiff(iRow,:);
    end
    
    subplot(1,3,iMod)
    yVals = nanmedian(cDiff);
    yVals = yVals(~isnan(yVals));
    xVals = xVals(~isnan(yVals));
    
    for iCol = 1:size(cDiff,2)
        data = cDiff(:,iCol);
        data = data(~isnan(data));
        nBoot = 1000;
        
        bootMedians = zeros(nBoot,1);
        for i = 1:nBoot
            resampleData = datasample(data,length(data));
            bootMedians(i) = median(resampleData);
        end
        ci = prctile(bootMedians,[2.5 97.5]);
        med = nanmedian(data);

        plot([iCol*5 iCol*5],[ci(1) ci(2)],'-','linewidth',1,'color',[0 0.4470 0.7410])
        hold on
%         plot(iMod,med,'.','markersize',20,'color',[0 0.4470 0.7410])
    end
    plot(xVals,yVals,'k.','markersize',15)
    hold on
    plot(xVals,yVals,'k-','markersize',15)

        hold on
    xlabel('Number of cells')
    ylabel('Rotation difference')
%     yMax = round(nanmax(rotMean))+0.5;
    ylim([-1 58])
    xlim([0 100])
    
end

%%
load distanceFinal_ed6b

figure('position',[-1778,528,1311,262]);
for iMod = 1:3
    distanceAll = nan(size(dFinal,1),25);
    for iRow = 1:size(dFinal,1)
        temp = dFinal{iRow,iMod};
        distanceAll(iRow,1:size(temp,1)) = temp;
    end
    
    xVals = 5:5:size(distanceAll,2)*5;
    for iRow = 1:size(distanceAll,1)
        yVals = distanceAll(iRow,:);
    end

    subplot(1,3,iMod)
        xVals = 5:5:size(distanceAll,2)*5;
    %     yVals = nanmean(distanceAll);
    
    yVals = nanmedian(distanceAll);
    yVals = yVals(~isnan(yVals));
    xVals = xVals(~isnan(yVals));
    
    for iCol = 1:size(distanceAll,2)
        data = distanceAll(:,iCol);
        data = data(~isnan(data));
        nBoot = 1000;
        
        bootMedians = zeros(nBoot,1);
        for i = 1:nBoot
            resampleData = datasample(data,length(data));
            bootMedians(i) = median(resampleData);
        end
        ci = prctile(bootMedians,[2.5 97.5]);
        med = nanmedian(data);

        plot([iCol*5 iCol*5],[ci(1) ci(2)],'-','linewidth',1,'color',[0 0.4470 0.7410])
        hold on
%         plot(iMod,med,'.','markersize',20,'color',[0 0.4470 0.7410])
    end
    plot(xVals,yVals,'k.','markersize',15)
    hold on
    plot(xVals,yVals,'k-','markersize',15)
    
    xlabel('Number of cells')
    ylabel('Distance')
    %     yMax = round(nanmax(rotMean))+0.5;
    ylim([-1 20])
    xlim([0 100])
    
end