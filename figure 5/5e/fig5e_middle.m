%%
load rotCCstore_5e
figure;

rotCC = rotCCstore(:,:,1);
maxCorr(1,1) = nanmax(nanmax(rotCC));
subplot(131)
rotCC(isnan(rotCC)) = eps;
r = general.smooth(rotCC,4);
colorMapBRK(r);
maxVal = nanmax(nanmax(r));
minVal = nanmin(nanmin(r));
caxis([minVal maxVal])
xline(size(rotCC,2)/2)
yline(size(rotCC,2)/2)

rotCC = rotCCstore(:,:,2);
maxCorr(1,2) = nanmax(nanmax(rotCC));
subplot(132)
rotCC(isnan(rotCC)) = eps;
r = general.smooth(rotCC,4);
colorMapBRK(r);
caxis([minVal maxVal])
xline(size(rotCC,2)/2)
yline(size(rotCC,2)/2)

rotCC = rotCCstore(:,:,3);
maxCorr(1,3) = nanmax(nanmax(rotCC))
subplot(133)
r = general.smooth(rotCC,4);
colorMapBRK(r);
caxis([minVal maxVal])
xline(size(rotCC,2)/2)
yline(size(rotCC,2)/2)
