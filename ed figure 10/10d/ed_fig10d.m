load fig3de_sessionPairs
load expList_HP_3de
groupIDs = [1 3 2];

load colorMat

%%
load pfShift_ed10d

figure('position',[-1719,38,1663,297]);
for iGrp = 1:3
    subplot(1,4,iGrp)
    inds = find(idx == groupIDs(iGrp));
    
    maxStore = nan(size(inds,1),1);
    for i = 1:size(inds,1)
        if ~isempty(expStore{inds(i),4})
            distStore = expStore{inds(i),4};
            
            pd = fitdist(distStore,'kernel');
            x = 0:56;
            y = pdf(pd,x);
            plot(x,y,'color',colorMat(iGrp,:))
            hold on
            
            [maxVal,maxInd] = nanmax(y);
            xlim([0 56])
            ylim([0 0.1])
            
            maxStore(i,1) = maxInd;
        end
    end
    
    grid off
    title ''
    xlabel('Place field shift')
    ylabel('Frequency')
    axis square
    
    grpStore{1,iGrp} = maxStore;
end

subplot(1,4,4)
maxStore = nan(size(expStoreStab,1),1);

for i = 1:size(expStoreStab,1)
    distStore = expStoreStab{i,4};
    
    if ~isempty(distStore)
        pd = fitdist(distStore,'kernel');
        x = 0:56;
        y = pdf(pd,x);
        plot(x,y,'color',[0.5 0.5 0.5])
        hold on
        
        [maxVal,maxInd] = nanmax(y);

        xlim([0 56])
        ylim([0 0.22])
        
        maxStore(i,1) = maxInd;
    end
end

grid off
title ''
xlabel('Place field shift')
ylabel('Frequency')
axis square

grpStore{1,4} = maxStore;

%%
figure('position',[-863,39,277,296]);
plotSpread({grpStore{1,1} grpStore{1,2} grpStore{1,3} grpStore{1,4}})
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','o','markersize',6,'linewidth',1.5)
box on
ylabel('Place field shift (cm)')
ylim([0 32])

for iGrp = 1:4
    grpMed(1,iGrp) = nanmedian(grpStore{1,iGrp});
end
plot(1:4,grpMed,'r+','markersize',15,'linewidth',2)
set(gca,'ytick',0:10:50)
xlabel('Group')

