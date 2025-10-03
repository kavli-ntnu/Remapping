% load rotChange_ed11h

load colorVec
figure;
scatter(nanmean(rotChange,2),meanStore,35,colorVec)
xlim([-185 185])
set(gca,'xtick',-180:60:180)
box on
ylim([-0.02 0.72])
xlabel('Rotation (deg)')
ylabel('Spatial correlation')

%%
groupIDs = [2 3 1];
rotChange = nanmean(abs(rotChange),2);
grp1rot = rotChange(idx == groupIDs(1));
grp1sc = meanStore(idx == groupIDs(1));

figure;
scatter(grp1rot,grp1sc,35)
xlim([-5 185])
ylim([-0.02 0.72])
xlabel('Rotation (deg)')
ylabel('Spatial correlation')

xVals = grp1rot;
yVals = grp1sc;
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);

[r,p] = corr(xVals,yVals)
