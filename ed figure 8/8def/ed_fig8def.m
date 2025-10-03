% CA1
load storeCA1_ed8
placeCorrStore = placeCorrStoreCA1;
rateDiffStore = rateDiffStoreCA1;

figure('position',[-1552,147,1301,451]);

medianStore = nan(size(placeCorrStore,1),4);
medianRateStore = nan(size(placeCorrStore,1),4);

remapMat = [];
stabMat = [];

remapRate = [];
stabRate = [];

for iRow = 1:size(placeCorrStore,1) 
    if ~isempty(placeCorrStore{iRow,1})
        temp = placeCorrStore{iRow,1};
        medianStore(iRow,1:size(temp,2)) = nanmedian(temp);
        
        tempRate = rateDiffStore{iRow,1};
        medianRateStore(iRow,1:size(tempRate,2)) = nanmedian(tempRate);
        
        if size(temp,2) == 4
            subplot(121)
            C = cdfplot(temp(:,3));
            set(C,'color',[153/255 0 0])
            hold on
            
            remapMat = [remapMat; temp(:,3)];
            stabMat = [stabMat; temp(:,1)];
            
            subplot(122)
            C = cdfplot(tempRate(:,3))
            set(C,'color',[153/255 0 0])
            hold on
            
            remapRate = [remapRate; tempRate(:,3)];
            stabRate = [stabRate; tempRate(:,1)];
            
        elseif size(temp,2) == 2
            subplot(121)
            C = cdfplot(temp(:,2));
            set(C,'color',[153/255 0 0])
            hold on
            
            remapMat = [remapMat; temp(:,2)];
            
            subplot(122)
            C = cdfplot(tempRate(:,1));
            set(C,'color',[153/255 0 0])
            hold on
            
            remapRate = [remapRate; tempRate(:,1)];
        end
    end
end

subplot(121)
box on
title ''
grid off
xlabel('Spatial correlation')
ylabel('Proportion')
set(gca,'xtick',-0.5:0.5:1)
set(gca,'ytick',0:0.2:1)
xlim([-0.5 1])
title('CA1')

% hold on
% cdfplot(remapMat)
% cdfplot(stabMat)

subplot(122)
box on
title ''
grid off
xlabel('Rate difference')
ylabel('Proportion')
set(gca,'xtick',0:0.2:1)
set(gca,'ytick',0:0.2:1)
title('CA1')

remapMatCA1 = remapMat;
stabMatCA1 = stabMat;
remapRateCA1 = remapRate;
stabRateCA1 = stabRate;

%% CA3
load storeCA3_ed8
placeCorrStore = placeCorrStoreCA3;
rateDiffStore = rateDiffStoreCA3;

figure('position',[-1552,147,1301,451]);

medianStore = nan(size(placeCorrStore,1),4);
medianRateStore = nan(size(placeCorrStore,1),4);

remapMat = [];
stabMat = [];

remapRate = [];
stabRate = [];

for iRow = 1:size(placeCorrStore,1) 
    if ~isempty(placeCorrStore{iRow,1})
        temp = placeCorrStore{iRow,1};
        medianStore(iRow,1:size(temp,2)) = nanmedian(temp);
        
        tempRate = rateDiffStore{iRow,1};
        medianRateStore(iRow,1:size(tempRate,2)) = nanmedian(tempRate);
        
        if size(temp,2) == 4
            subplot(121)
            C = cdfplot(temp(:,3));
            set(C,'color',[0 102/255 204/255])
            hold on
            
            remapMat = [remapMat; temp(:,3)];
            stabMat = [stabMat; temp(:,1)];
            
            subplot(122)
            C = cdfplot(tempRate(:,3));
            set(C,'color',[0 102/255 204/255])
            hold on
            
            remapRate = [remapRate; tempRate(:,3)];
            stabRate = [stabRate; tempRate(:,1)];
            
        elseif size(temp,2) == 2
            subplot(121)
            C = cdfplot(temp(:,2));
            set(C,'color',[0 102/255 204/255])
            hold on
            
            remapMat = [remapMat; temp(:,2)];
            
            subplot(122)
            C = cdfplot(tempRate(:,1));
            set(C,'color',[0 102/255 204/255])
            hold on
            
            remapRate = [remapRate; tempRate(:,1)];
        end
    end
end

subplot(121)
box on
title ''
grid off
xlabel('Spatial correlation')
ylabel('Proportion')
set(gca,'xtick',-0.5:0.5:1)
set(gca,'ytick',0:0.2:1)
xlim([-0.5 1])
title('CA3')

% hold on
% cdfplot(remapMat)
% cdfplot(stabMat)

subplot(122)
box on
title ''
grid off
xlabel('Rate difference')
ylabel('Proportion')
set(gca,'xtick',0:0.2:1)
set(gca,'ytick',0:0.2:1)
title('CA3')

remapMatCA3 = remapMat;
stabMatCA3 = stabMat;
remapRateCA3 = remapRate;
stabRateCA3 = stabRate;

%%
figure;
h1 = histfit(remapMatCA1,100,'kernel');
x1 = h1(2).XData;
y1 = h1(2).YData;
counts1 = histcounts(remapMatCA1,100);
close(gcf)

figure;
h2 = histfit(stabMatCA1,100,'kernel');
x2 = h2(2).XData;
y2 = h2(2).YData;
counts2 = histcounts(stabMatCA1,100);
close(gcf)

figure;
h3 = histfit(remapMatCA3,100,'kernel');
x3 = h3(2).XData;
y3 = h3(2).YData;
counts3 = histcounts(remapMatCA3,100);
close(gcf)

figure;
h4 = histfit(stabMatCA3,100,'kernel');
x4 = h4(2).XData;
y4 = h4(2).YData;
counts4 = histcounts(stabMatCA3,100);
close(gcf)

%% plot them together!
figure;
plot(x1,y1/sum(counts1),'-','color',[153/255 0 0])
hold on
plot(x2,y2/sum(counts2),'-','color',[1 153/255 153/255])
plot(x3,y3/sum(counts3),'-','color',[0 102/255 204/255])
plot(x4,y4/sum(counts4),'-','color',[153/255 204/255 1])
xlabel('Spatial correlation')
ylabel('Frequency')

xlim([-0.5 1])
set(gca,'ytick',0:0.01:0.04)

[max1,mind1] = nanmax(y1/sum(counts1));
[max2,mind2] = nanmax(y2/sum(counts2));
[max3,mind3] = nanmax(y3/sum(counts3));
[max4,mind4] = nanmax(y4/sum(counts4));

%% REPEAT FOR RATE
figure;
h1 = histfit(remapRateCA1,100,'kernel');
x1 = h1(2).XData;
y1 = h1(2).YData;
counts1 = histcounts(remapRateCA1,100);
close(gcf)

figure;
h2 = histfit(stabRateCA1,100,'kernel');
x2 = h2(2).XData;
y2 = h2(2).YData;
counts2 = histcounts(stabRateCA1,100);
close(gcf)

figure;
h3 = histfit(remapRateCA3,100,'kernel');
x3 = h3(2).XData;
y3 = h3(2).YData;
counts3 = histcounts(remapRateCA3,100);
close(gcf)

figure;
h4 = histfit(stabRateCA3,100,'kernel');
x4 = h4(2).XData;
y4 = h4(2).YData;
counts4 = histcounts(stabRateCA3,100);
close(gcf)

%% plot them together!
figure;
plot(x1,y1/sum(counts1),'-','color',[153/255 0 0])
hold on
plot(x2,y2/sum(counts2),'-','color',[1 153/255 153/255])
plot(x3,y3/sum(counts3),'-','color',[0 102/255 204/255])
plot(x4,y4/sum(counts4),'-','color',[153/255 204/255 1])
xlabel('Rate difference')
ylabel('Frequency')

xlim([-0.1 1])
ylim([0 0.027])
set(gca,'xtick',-0.2:0.2:1,'ytick',0:0.01:0.04)

[max1,mind1] = nanmax(y1/sum(counts1));
[max2,mind2] = nanmax(y2/sum(counts2));
[max3,mind3] = nanmax(y3/sum(counts3));
[max4,mind4] = nanmax(y4/sum(counts4));
