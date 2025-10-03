load grpRDstore_ed10g
load pvCorr_ed10e

load colorMat
colorMat(4,:) = [0.5 0.5 0.5];

figure('position',[-1719,38,1663,297]);
for iGrp = 1:4
    subplot(1,4,iGrp)
   
    currentGrp = grpRDstore{1,iGrp};
    for i = 1:size(currentGrp,1)
        if ~isempty(currentGrp{i,1})
            C = cdfplot(currentGrp{i,1});
            hold on
            set(C,'color',colorMat(iGrp,:))
        end
    end
    grid off
    title ''
    xlabel('Rate difference')
    ylabel('Frequency')
end

%%
figure;
for iGrp = 1:4
    grpVals = cell2mat(grpRDstore{1,iGrp});
    
    [f,x,flo,fup] = ecdf(grpVals);
    shadedplot(x',flo',fup',colorMat(iGrp,:),colorMat(iGrp,:));
    C = cdfplot(grpVals);
    set(C,'color','k','linewidth',0.5)
    hold on
    
    grid off
    title ''
    xlabel('PV correlation')
    ylabel('Frequency')
end

%%
figure;
m1_col = [];
m2_col = [];
e2_col = [];
for iGrp = 1:3
    
    m1 = grpMinDist{1,iGrp};
    m2 = grpRDstore{1,iGrp};
    
    m2_meanStore = nan(size(m2,1),1);
    e2_meanStore = nan(size(m2,1),1);
    for iRow = 1:size(m2,1)
        m2_meanStore(iRow,1) = nanmedian(m2{iRow,1});
        e2_errStore(iRow,1) = nanstd(m2{iRow,1}) ./ sqrt(sum(~isnan(m2{iRow,1})));
    end
    
    m1_col = [m1_col; m1];
    m2_col = [m2_col; m2_meanStore];
    e2_col = [e2_col; e2_errStore];
    
    scatter(m1,m2_meanStore,35,colorMat(iGrp,:),'linewidth',1)
    hold on
end

for iPt = 1:size(m2_col,1)
    errorbar(m1_col(iPt,1),m2_col(iPt,1),e2_col(iPt,1),'capsize',0,'color',[0.8 0.8 0.8],'linewidth',1.2)
    hold on
end
    
xlim([-0.02 0.72])
ylim([0 0.45])
set(gca,'xtick',0:0.1:0.8)
xlabel('Minimum distance')
ylabel('Mean rate change')
box on

xVals = m1_col;
yVals = m2_col;
xLogic = ~isnan(xVals);
yLogic = ~isnan(yVals);
xVals = xVals(xLogic & yLogic);
yVals = yVals(xLogic & yLogic);
[r,p] = corr(xVals,yVals)
