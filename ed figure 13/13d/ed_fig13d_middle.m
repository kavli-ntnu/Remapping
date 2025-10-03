load rotCCstore_rotShift_ed13d

figure;
subplot(121)
rotCC = rotCCstore(:,:,1);
rotCC(isnan(rotCC)) = eps;
r = general.smooth(rotCC,4);
colorMapBRK(r);
xline(size(rotCC,2)/2)
yline(size(rotCC,2)/2)
c = caxis;
title('Sim 5617')

subplot(122)
rotCC = rotCCstore(:,:,2);
rotCC(isnan(rotCC)) = eps;
r = general.smooth(rotCC,4);
colorMapBRK(r);
caxis(c)
xline(size(rotCC,2)/2)
yline(size(rotCC,2)/2)
title('Sim 7438')
