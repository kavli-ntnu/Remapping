load mxAA_mxAB_ed8b

%%
mxAAerr = nanstd(mxAA) ./ sqrt(sum(~isnan(mxAA)));
mxABerr = nanstd(mxAB) ./ sqrt(sum(~isnan(mxAB)));

figure;
[ha hb hc] = shadedplot(1:120,nanmean(mxAB)+mxAAerr,nanmean(mxAB)-mxAAerr,[128/225 128/225 128/225],[0 0 0]);
hold on
plot(1:120,nanmean(mxAB),'k-','linewidth',2)

[ha hb hc] = shadedplot(1:120,nanmean(mxAA)+mxABerr,nanmean(mxAA)-mxABerr,[224/255 224/255 224/255],[0 0 0]);
hold on
plot(1:120,nanmean(mxAA),'k-','linewidth',2)

grid off
ylabel('PV correlation')
set(gca,'xtick',0:20:120,'xticklabel',0:60:360)
xlabel('Rotation')
ylim([0 0.7])
set(gca,'ytick',0:0.2:1)

%% at optimal rotation
maxAA = nanmax(mxAA,[],2);
maxAB = nanmax(mxAB,[],2);
   
figure('position',[1285,67,198,420]);
plotSpread({maxAA(:,1) maxAB(:,1)})
set(findobj(gca,'type','line'),'color','k','marker','o','markersize',6,'linewidth',1.5)
hold on
plot(1:2,[nanmean(maxAA(:,1)) nanmean(maxAB(:,1))],'r+','markersize',15,'linewidth',2)
set(gca,'xtick',1:2)
xlim([0.5 2.5])
ylim([0 0.7])
set(gca,'ytick',0:0.2:1)
ylabel('PV correlation')
