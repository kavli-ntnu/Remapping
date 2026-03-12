load fig3de_sessionPairs
load expList_HP_3de
groupIDs = [1 3 2];

load colorMat
load mxStore_ed10e
grpStore = cell(1);

figure('position',[-1714,23,1614,282]);
for iGrp = 1:size(groupIDs,2)
    subplot(1,4,iGrp)
    inds = find(idx == groupIDs(iGrp));
    mxTemp = mxStoreAll(inds,:);
    grpStore{1,iGrp} = mxTemp;
    
    for iRow = 1:size(mxTemp,1)
        if inds(iRow) == 23
            m = mxTemp(iRow,:);
            [maxVal,ind] = nanmax(m);
            m = circshift(m,-ind);
            plot(1:120,m,'color',colorMat(iGrp,:))
        else
            plot(1:120,mxTemp(iRow,:),'color',colorMat(iGrp,:))
            hold on
        end
    end
    ylim([0 0.75])
    set(gca,'ytick',0:0.25:1,'xtick',0:20:120,'xticklabels',0:60:360)
    ylabel('PV correlation')
    xlabel('Rotation')
end

load mxStable_ed10e
grpStore{1,4} = mxStable;
subplot(1,4,4)

for iRow = 1:size(mxStable,1)
    plot(1:120,mxStable(iRow,:),'color',[0.5 0.5 0.5])
    hold on
end
ylim([0 0.75])
set(gca,'ytick',0:0.25:1,'xtick',0:20:120,'xticklabels',0:60:360)
ylabel('PV correlation')
xlabel('Rotation')

%%
mxGlobal = grpStore{1,1};
mxPartial = grpStore{1,2};
mxWeak = grpStore{1,3};
mxCtrl = grpStore{1,4};

mxGlobalErr = nanstd(grpStore{1,1}) ./ sqrt(sum(~isnan(grpStore{1,1})));
mxPartialErr = nanstd(grpStore{1,2}) ./ sqrt(sum(~isnan(grpStore{1,2})));
mxWeakErr = nanstd(grpStore{1,3}) ./ sqrt(sum(~isnan(grpStore{1,3})));
mxStableErr = nanstd(grpStore{1,4}) ./ sqrt(sum(~isnan(grpStore{1,4})));

figure;
[ha hb hc] = shadedplot(1:120,nanmean(mxGlobal)+mxGlobalErr,nanmean(mxGlobal)-mxGlobalErr,[128/225 128/225 128/225],[0 0 0]);
hold on
plot(1:120,nanmean(mxGlobal),'k-','linewidth',2)

[ha hb hc] = shadedplot(1:120,nanmean(mxPartial)+mxPartialErr,nanmean(mxPartial)-mxPartialErr,[180/255 180/255 180/255],[0 0 0]);
hold on
plot(1:120,nanmean(mxPartial),'k-','linewidth',2)

[ha hb hc] = shadedplot(1:120,nanmean(mxWeak)+mxWeakErr,nanmean(mxWeak)-mxWeakErr,[224/255 224/255 224/255],[0 0 0]);
hold on
plot(1:120,nanmean(mxWeak),'k-','linewidth',2)

[ha hb hc] = shadedplot(1:120,nanmean(mxStable)+mxStableErr,nanmean(mxStable)-mxStableErr,[240/255 240/255 240/255],[0 0 0]);
hold on
plot(1:120,nanmean(mxStable),'k-','linewidth',2)

grid off
ylabel('PV correlation')
set(gca,'xtick',0:20:120,'xticklabel',0:60:360)
xlabel('Rotation')
ylim([0 0.75])
set(gca,'ytick',0:0.25:1)
