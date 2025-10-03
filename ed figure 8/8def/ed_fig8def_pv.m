
load pvStore_ed8_def

%%
figure;
hold on
for iRow = 1:size(pvStoreAB_CA1,1)
    C = cdfplot(pvStoreAB_CA1{iRow,1});
    set(C,'color',[153/255 0 0])
end

title ''
grid off
xlabel('PV correlation')
ylabel('Proportion')
xlim([-0.5 1])
set(gca,'ytick',0:0.2:1)

%%
figure;
hold on
for iRow = 1:size(pvStoreAB_CA3,1)
    C = cdfplot(pvStoreAB_CA3{iRow,1});
    set(C,'color',[0 102/255 204/255])
end

title ''
grid off
xlabel('PV correlation')
ylabel('Proportion')
xlim([-0.5 1])
set(gca,'ytick',0:0.2:1)

%%
CA1_AB = cell2mat(pvStoreAB_CA1);
CA1_AA = cell2mat(pvStoreAA_CA1);

CA3_AB = cell2mat(pvStoreAB_CA3);
CA3_AA = cell2mat(pvStoreAA_CA3);

%%
figure;
h1 = histfit(CA1_AB,50,'kernel');
x1 = h1(2).XData;
y1 = h1(2).YData;
counts1 = histcounts(CA1_AB,50);
close(gcf)

figure;
h2 = histfit(CA1_AA,50,'kernel');
x2 = h2(2).XData;
y2 = h2(2).YData;
counts2 = histcounts(CA1_AA,50);
close(gcf)

figure;
h3 = histfit(CA3_AB,50,'kernel');
x3 = h3(2).XData;
y3 = h3(2).YData;
counts3 = histcounts(CA3_AB,50);
close(gcf)

figure;
h4 = histfit(CA3_AA,50,'kernel');
x4 = h4(2).XData;
y4 = h4(2).YData;
counts4 = histcounts(CA3_AA,50);
close(gcf)

%% plot them together!
figure;
plot(x1,y1/sum(counts1),'-','color',[153/255 0 0])
hold on
plot(x2,y2/sum(counts2),'-','color',[1 153/255 153/255])
plot(x3,y3/sum(counts3),'-','color',[0 102/255 204/255])
plot(x4,y4/sum(counts4),'-','color',[153/255 204/255 1])
xlabel('PV correlation')
ylabel('Frequency')

xlim([-0.5 1.1])
set(gca,'ytick',0:0.02:0.08)

[max1,mind1] = nanmax(y1/sum(counts1));
[max2,mind2] = nanmax(y2/sum(counts2));
[max3,mind3] = nanmax(y3/sum(counts3));
[max4,mind4] = nanmax(y4/sum(counts4));

