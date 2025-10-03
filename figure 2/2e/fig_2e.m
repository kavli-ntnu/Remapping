load realignmentPV_29730

%%
figure;
subplot(121)
colorMapBRK(rotCCstore(:,:,mindStore(1),1))
caxis([-inf maxCorr(1,1)])
hold on
xline(size(rotCCstore,1)/2)
yline(size(rotCCstore,1)/2)

subplot(122)
colorMapBRK(rotCCstore(:,:,mindStore(2),2))
caxis([-inf maxCorr(1,1)])
hold on
xline(size(rotCCstore,1)/2)
yline(size(rotCCstore,1)/2)

maxCorr
rotAngle