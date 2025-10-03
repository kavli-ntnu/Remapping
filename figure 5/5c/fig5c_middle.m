load shiftPV

%%
simNum = 1382;
rotCC = rotCCstore(:,:,simNum);

figure;
subplot(131)
r = general.smooth(rotCC,4);
colorMapBRK(r);
minVal = nanmin(nanmin(r));
maxVal = nanmax(nanmax(r));
caxis([minVal maxVal])
xline(size(rotCC,2)/2)
yline(size(rotCC,2)/2)
title(sprintf('%.2f',maxCorr(simNum)))

simNum = 4858;
rotCC = rotCCstore(:,:,simNum);

subplot(132)
r = general.smooth(rotCC,4);
colorMapBRK(r);
caxis([minVal maxVal])
xline(size(rotCC,2)/2)
yline(size(rotCC,2)/2)
title(sprintf('%.2f',maxCorr(simNum)))


simNum = 7347;
rotCC = rotCCstore(:,:,simNum);

subplot(133)
r = general.smooth(rotCC,4);
colorMapBRK(r);
caxis([minVal maxVal])
xline(size(rotCC,2)/2)
yline(size(rotCC,2)/2)
title(sprintf('%.2f',maxCorr(simNum)))
