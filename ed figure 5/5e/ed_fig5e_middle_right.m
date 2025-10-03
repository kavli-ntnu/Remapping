load numPerMod_ed5e
load distStore_ed5e

numPair(:,1) = numPerMod(:,1) + numPerMod(:,2);
numPair(:,2) = numPerMod(:,1) + numPerMod(:,3);
numPair(:,3) = numPerMod(:,2) + numPerMod(:,3);

n = reshape(numPair,[],1);
d = reshape(distStore,[],1);

figure;
scatter(n,d)
box on
xlabel('Number of grid cells')
xlim([0 275])
ylim([0 90])
set(gca,'ytick',0:20:100)
ylabel('Distance (cm)')

xVals = n;
yVals = d;
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);

[r,p] = corr(xVals,yVals)

%%
load maxCC_ed5e

pairCC = [];
pairCC(:,1) = nanmean(maxCC(:,1:2),2);
pairCC(:,2) = nanmean(maxCC(:,[1 3]),2);
pairCC(:,3) = nanmean(maxCC(:,[2 3]),2);

minCC = [];
minCC(:,1) = nanmin(maxCC(:,1:2),[],2);
minCC(:,2) = nanmin(maxCC(:,[1 3]),[],2);
minCC(:,3) = nanmin(maxCC(:,[2 3]),[],2);

minCC = reshape(minCC,[],1);
pairCC = reshape(pairCC,[],1);
pairDist = reshape(distStore,[],1);

figure;
scatter(minCC,pairDist)
box on
xlabel('Min PV correlation')
xlim([0.08 0.38])
ylim([0 90])
set(gca,'ytick',0:20:100)
ylabel('Distance (cm)')

xVals = minCC;
yVals = pairDist;
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);

[r,p] = corr(xVals,yVals)

