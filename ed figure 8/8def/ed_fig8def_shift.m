load pfShift_ed8_def

%%
figure('position',[-1552,147,1301,451]);

subplot(121)
distStoreCA1 = [];
for iExp = 1:size(expStoreCA1_remap,1)
    C = cdfplot(expStoreCA1_remap{iExp,1});
    set(C,'color',[153/255 0 0])
    hold on    
    distStoreCA1 = [distStoreCA1; expStoreCA1_remap{iExp,1}];
end

box on
title ''
grid off
xlabel('Place field shift (cm)')
ylabel('Proportion')
set(gca,'xtick',0:20:60)
set(gca,'ytick',0:0.2:1)
title('CA1')

subplot(122)
distStoreCA3 = [];
for iExp = 1:size(expStoreCA3_remap,1)
    C = cdfplot(expStoreCA3_remap{iExp,1});
    set(C,'color',[0 102/255 204/255])
    hold on
    distStoreCA3 = [distStoreCA3; expStoreCA3_remap{iExp,1}];
end

box on
title ''
grid off
xlabel('Place field shift (cm)')
ylabel('Proportion')
set(gca,'xtick',0:20:60)
set(gca,'ytick',0:0.2:1)
title('CA3')

%%
distStoreCA1_stab = [];
for iExp = 1:size(expStoreCA1_stab,1)
    distStoreCA1_stab = [distStoreCA1_stab; expStoreCA1_stab{iExp,1}];
end
distStoreCA3_stab = [];
for iExp = 1:size(expStoreCA3_stab,1)
    distStoreCA3_stab = [distStoreCA3_stab; expStoreCA3_stab{iExp,1}];
end

figure;
h1 = histfit(distStoreCA1,100,'kernel');
x1 = h1(2).XData;
y1 = h1(2).YData;
counts1 = histcounts(distStoreCA1,100);
close(gcf)

figure;
h2 = histfit(distStoreCA1_stab,100,'kernel');
x2 = h2(2).XData;
y2 = h2(2).YData;
counts2 = histcounts(distStoreCA1_stab,100);
close(gcf)

figure;
h3 = histfit(distStoreCA3,100,'kernel');
x3 = h3(2).XData;
y3 = h3(2).YData;
counts3 = histcounts(distStoreCA3,100);
close(gcf)

figure;
h4 = histfit(distStoreCA3_stab,100,'kernel');
x4 = h4(2).XData;
y4 = h4(2).YData;
counts4 = histcounts(distStoreCA3_stab,100);
close(gcf)

%% plot them together!
% subplot(122)

figure;
plot(x1,y1/sum(counts1),'-','color',[153/255 0 0])
hold on
plot(x2,y2/sum(counts2),'-','color',[1 153/255 153/255])
plot(x3,y3/sum(counts3),'-','color',[0 102/255 204/255])
plot(x4,y4/sum(counts4),'-','color',[153/255 204/255 1])
xlabel('Place field shift (cm)')
ylabel('Frequency')

xlim([-5 60])
set(gca,'xtick',0:20:60,'ytick',0:0.02:0.06)

[max1,mind1] = nanmax(y1/sum(counts1));
[max2,mind2] = nanmax(y2/sum(counts2));
[max3,mind3] = nanmax(y3/sum(counts3));
[max4,mind4] = nanmax(y4/sum(counts4));


