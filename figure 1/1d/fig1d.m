load expList_26821
szmap = 40;

load rotCCstore_AA_26821
numMods = size(rotCCstoreAA,4);

load pkStore_26821
pkStoreAB = pkStore{1,1};

load colorVec_ed1c
colorVec = colorVec / 255;
load colorMap

%
ptStore = cell(1);
ptShift = cell(1);
 
figure('position',[919,-49,1920,955]);
counter = 1;
numCols = 4;

for iMod = 1:numMods
    subplot(numMods,numCols,counter)
    colorMapBRK(rotCCstoreAA(:,:,mindStore(iMod),iMod));
    title(sprintf('%.2f',maxCorr(1,iMod)))
    
    lineVal = size(rotCCstoreAA,1) / 2;
    xline(lineVal)
    yline(lineVal)

    acTemp = rotCCstoreAA(:,:,mindStore(iMod),iMod);
    acTemp(isnan(acTemp)) = 0;
    [gScore,gStats,~,COM] = analyses.gridnessScore(acTemp,'threshold',0.05);
    comStore(:,:,1) = COM;

    acTemp(isnan(acTemp)) = 0;
    I = mat2gray(acTemp,[0 1]);
    strelNum = 1;
    SE = strel('sphere',strelNum);
    J = imerode(I,SE);
    IM = imreconstruct(J,acTemp);
    regmax = imregionalmax(IM);
    doubleMap = double(regmax);
    [f1,f2] = analyses.placefield(doubleMap,'minBins',5);

    ptsA1 = nan(size(f2,2),2);
    if ~isempty(f2)
        for iField = 1:size(f2,2)
            ptsA1(iField,1) = f2(iField).peakX; 
            ptsA1(iField,2) = f2(iField).peakY; 
        end
    end
    
    ptStore{iMod,1} = ptsA1;

    counter = counter + numCols;
end

%
load rotCCstore_AB_26821
rotCCstoreAB = rotCCstore;
counter = 2;

for iMod = 1:numMods
    subplot(numMods,numCols,counter)
    colorMapBRK(rotCCstoreAB(:,:,mindStore(iMod),iMod));
    title(sprintf('%.2f',maxCorr(1,iMod)))
    
    lineVal = size(rotCCstoreAB,1) / 2;
    xline(lineVal)
    yline(lineVal)
    
    counter = counter + numCols;
end

%
counter = 3;
numMods = 3;
for iMod = 1:numMods

    subplot(numMods,numCols,counter)
    ptTemp = ptStore{iMod,1};
    plot(ptTemp(:,1),ptTemp(:,2),'.','markersize',15,'color',colorVec(iMod,:,2))
    hold on
      
    currentPk = pkStoreAB(iMod,:);
    
    origin = szmap;
    changePk = nan(1,2);
    changePk(2) = origin-currentPk(2);
    changePk(1) = origin-currentPk(1);
    
    shiftTemp = ptTemp - changePk;
    
    hold on
    plot(shiftTemp(:,1),shiftTemp(:,2),'.','markersize',15,'color',colorVec(iMod,:,1))
    
    plot([origin pkStoreAB(iMod,1)],[origin pkStoreAB(iMod,2)],'k-')
    
    set(gca,'ydir','rev') % reverse to match field map
    axis square
    box on
    set(gca,'xticklabels','','yticklabels','')
    
    counter = counter + numCols;
    
    ptShift{iMod,1} = shiftTemp;
    
    inc = 5;
    xlim([0+inc size(acTemp,1)-inc])
    ylim([0+inc size(acTemp,1)-inc])

    lineVal = size(rotCCstoreAB,1) / 2;
    xline(lineVal)
    yline(lineVal)
end

