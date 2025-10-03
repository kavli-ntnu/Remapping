load expList_26718
load rotCCstore_26718

numMods = size(rotAngle,2);
figure;
for iMod = 1:numMods
    subplot(1,numMods,iMod)
    colorMapBRK(rotCCstore(:,:,mindStore(iMod),iMod));
    xline(40)
    yline(40)
    title(sprintf('%d deg',rotAngle(1,iMod)))
    text(0,85,sprintf('%.2f',maxCorr(1,iMod)))
end

%%
load run1_26718

figure;
for iMod = 1:3
    rotCCstore = rotCCstore_run1{1,iMod}.rotCCstore;
    mindStore = rotCCstore_run1{1,iMod}.mindStore;
    rotAngle = rotCCstore_run1{1,iMod}.rotAngle;
    maxCorr = rotCCstore_run1{1,iMod}.maxCorr;
    
    subplot(3,3,iMod) 
    colorMapBRK(rotCCstore(:,:,mindStore(iMod),iMod));
    hold on
    xline(40)
    yline(40)
    text(0,-10,sprintf('%d deg',rotAngle(iMod)))
    text(60,-10,sprintf('%.2f',maxCorr(iMod)))
end

%%
load run2_26718

for iMod = 1:3
    rotCCstore = rotCCstore_run2{1,iMod}.rotCCstore;
    mindStore = rotCCstore_run2{1,iMod}.mindStore;
    rotAngle = rotCCstore_run2{1,iMod}.rotAngle;
    maxCorr = rotCCstore_run2{1,iMod}.maxCorr;
    
    subplot(3,3,iMod+3) 
    colorMapBRK(rotCCstore(:,:,mindStore(iMod),iMod));
    hold on
    xline(40)
    yline(40)
    text(0,-10,sprintf('%d deg',rotAngle(iMod)))
    text(60,-10,sprintf('%.2f',maxCorr(iMod)))
end

%%
load run3_26718

for iMod = 1:3
    rotCCstore = rotCCstore_run3{1,iMod}.rotCCstore;
    mindStore = rotCCstore_run3{1,iMod}.mindStore;
    rotAngle = rotCCstore_run3{1,iMod}.rotAngle;
    maxCorr = rotCCstore_run3{1,iMod}.maxCorr;
    
    subplot(3,3,iMod+6) 
    colorMapBRK(rotCCstore(:,:,mindStore(iMod),iMod));
    hold on
    xline(40)
    yline(40)
    text(0,-10,sprintf('%d deg',rotAngle(iMod)))
    text(60,-10,sprintf('%.2f',maxCorr(iMod)))
end
