load pvCorr_ed10e

load colorMat
figure;

m1_col = [];
m2_col = [];
for iGrp = 1:3  
    m1 = grpMinDist{1,iGrp};    
    m2 = grpPVrealign{1,iGrp};

    scatter(m1,m2,35,colorMat(iGrp,:),'linewidth',1)
    hold on
 
    m1_col = [m1_col; m1];
    m2_col = [m2_col; m2];
end

ylim([-0.02 0.62])
xlim([-0.02 0.72])
set(gca,'xtick',0:0.1:1)
xlabel('Minimum distance')
ylabel('PV correlation')
box on

xVals = m1_col;
yVals = m2_col;
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);
[r,p] = corr(xVals,yVals)
