% 3C, LEFT
load expList_HP_3c
load stabCorrStore_3c

%%
figure('position',[-617,453,560,194]);
subplot(131)
plotSpread(expList_HP{1,8})
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','o','markersize',4,'linewidth',1)
hold on
plot(1,nanmedian(expList_HP{1,8}),'r+','markersize',15,'linewidth',2)
ylim([-0.4 1])
set(gca,'xtick','')
ylabel('Spatial correlation')
set(gca,'ytick',0:0.5:15)

subplot(132)
plotSpread(expList_HP{3,8})
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','o','markersize',4,'linewidth',1)
hold on
plot(1,nanmedian(expList_HP{3,8}),'r+','markersize',15,'linewidth',2)
ylim([-0.4 1])
set(gca,'xtick','')
set(gca,'ytick',0:0.5:1)

subplot(133)
plotSpread(expList_HP{2,8})
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','o','markersize',4,'linewidth',1)
hold on
plot(1,nanmedian(expList_HP{2,8}),'r+','markersize',15,'linewidth',2)
ylim([-0.4 1])
set(gca,'xtick','')
set(gca,'ytick',0:0.5:1)

%%
bw = 0.12;

figure('position',[-487,330,315,274])
subplot(311)
[f,xi,u] = ksdensity(expList_HP{1,8},-0.4:0.05:1,'bandwidth',bw)
[pks,loc,width,prom] = findpeaks(f,xi);
findpeaks(f,xi)
hold on

[f,xi,u] = ksdensity(stabCorrStore{1,1}(:,1),-0.4:0.05:1,'bandwidth',bw)
[pks,loc,width,prom] = findpeaks(f,xi);
findpeaks(f,xi)
grid off
xlim([-0.5 1])
box off

subplot(312)
[f,xi,u] = ksdensity(expList_HP{3,8},-0.4:0.05:1,'bandwidth',bw)
[pks,loc,width,prom] = findpeaks(f,xi);
findpeaks(f,xi)
hold on

[f,xi,u] = ksdensity(stabCorrStore{3,1}(:,1),-0.4:0.05:1,'bandwidth',bw)
[pks,loc,width,prom] = findpeaks(f,xi);
findpeaks(f,xi)
grid off
xlim([-0.5 1])
box off
ylabel('Frequency')

subplot(313)
[f,xi,u] = ksdensity(expList_HP{2,8},-0.4:0.05:1,'bandwidth',0.11)
[pks,loc,width,prom] = findpeaks(f,xi);
findpeaks(f,xi)
hold on

[f,xi,u] = ksdensity(stabCorrStore{2,1}(:,1),-0.4:0.05:1,'bandwidth',bw)
[pks,loc,width,prom] = findpeaks(f,xi);
findpeaks(f,xi)
grid off
xlim([-0.5 1])
box off
xlabel('Spatial correlation')

%%
% 3C, MIDDLE
load realignPV_HP_28367_0922
uiopen('rotCC_HP_28367_0922.fig',1)
caxis([-inf maxCorr])

uiopen('rotCC_HP_28721_0703.fig',1)
caxis([-inf maxCorr])

uiopen('rotCC_HP_28367_0923.fig',1)
caxis([-inf maxCorr])

%%
% 3C, RIGHT
load pvStore_28367_0923

figure('position',[-740,170,560,275]);
subplot(131)
plotSpread(pvStore)
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','.','markersize',6,'linewidth',1.5)
hold on
plot(1,nanmedian(pvStore),'r+','markersize',15,'linewidth',2)
set(gca,'xtick','','ytick',-0.2:0.2:1)
ylim([-0.2 0.9])
ylabel('PV correlation')

load pvStore_28721_0703

subplot(132)
plotSpread(pvStore)
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','.','markersize',6,'linewidth',1.5)
hold on
plot(1,nanmedian(pvStore),'r+','markersize',15,'linewidth',2)
set(gca,'xtick','','ytick',-0.2:0.2:1)
ylim([-0.2 0.9])

load pvStore_28367_0922

subplot(133)
plotSpread(pvStore)
set(findobj(gca,'type','line'),'color',[0.5 0.5 0.5],'marker','.','markersize',6,'linewidth',1.5)
hold on
plot(1,nanmedian(pvStore),'r+','markersize',15,'linewidth',2)
set(gca,'xtick','','ytick',-0.2:0.2:1)
ylim([-0.2 0.9])

%%
load pvStore_28367_0923

figure;
h1 = histfit(pvStore,100,'kernel');
x1 = h1(2).XData;
y1 = h1(2).YData;
counts1 = histcounts(pvStore,100);
close(gcf)

load pvStoreStab_28367_0922 % no A-A comparison here, use previous day

figure;
h2 = histfit(pvStoreStab,100,'kernel');
x2 = h2(2).XData;
y2 = h2(2).YData;
counts2 = histcounts(pvStoreStab,100);
close(gcf)

figure('position',[-502,325,137,299]);
plot(x1,y1/sum(counts1),'-')
hold on
plot(x2,y2/sum(counts2),'-')

xlim([-0.22 1.02])
camroll(-90)
ylim([0 0.035])
xlabel('PV correlation')

%%
load pvStore_28721_0703

figure;
h3 = histfit(pvStore,100,'kernel');
x3 = h3(2).XData;
y3 = h3(2).YData;
counts3 = histcounts(pvStore,100);
close(gcf)

load pvStoreStab_28721_0703

figure;
h4 = histfit(pvStoreStab,100,'kernel');
x4 = h4(2).XData;
y4 = h4(2).YData;
counts4 = histcounts(pvStoreStab,100);
close(gcf)

figure('position',[-502,325,137,299]);
plot(x3,y3/sum(counts3),'-')
hold on
plot(x4,y4/sum(counts4),'-')

xlim([-0.22 1.02])
camroll(-90)
ylim([0 0.035])
xlabel('PV correlation')

%%
load pvStore_28367_0922

figure;
h5 = histfit(pvStore,100,'kernel');
x5 = h5(2).XData;
y5 = h5(2).YData;
counts5 = histcounts(pvStore,100);
close(gcf)

load pvStoreStab_28367_0922

figure;
h6 = histfit(pvStoreStab,100,'kernel');
x6 = h6(2).XData;
y6 = h6(2).YData;
counts6 = histcounts(pvStoreStab,100);
close(gcf)

figure('position',[-502,325,137,299]);
plot(x5,y5/sum(counts5),'-')
hold on
plot(x6,y6/sum(counts6),'-')

xlim([-0.22 1.02])
camroll(-90)
ylim([0 0.035])
xlabel('PV correlation')

