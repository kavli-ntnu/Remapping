load pkDiff_rotDiff_ed6_defg

runString = {'d1','d2'};

for iRun = 1:2
    
    if iRun == 1 || iRun == 3
        figure('position',[-1042,15,991,853]);
        cnt = 1;
    end
    
    subplot(3,3,cnt)
    for iMod = 1:3
        data = pkDiffNormFinal{iRun,1}(:,iMod);
        data = data(~isnan(data));
        nBoot = 1000;
        
        bootMedians = zeros(nBoot,1);
        for i = 1:nBoot
            resampleData = datasample(data,length(data));
            bootMedians(i) = median(resampleData);
        end
        ci = prctile(bootMedians,[2.5 97.5]);
        med = nanmedian(data);

        plot([iMod iMod],[ci(1) ci(2)],'-','linewidth',1,'color',[0 0.4470 0.7410])
        hold on
        plot(iMod,med,'.','markersize',20,'color',[0 0.4470 0.7410])
    end
    
    temp = pkDiffNormFinal{iRun,1};
    if size(temp,2) > 3
        set(gca,'xtick',1:size(temp,2),'xticklabels',{'M1','M2','M3','M4'})  
    else
        set(gca,'xtick',1:size(temp,2),'xticklabels',{'M1','M2','M3'})
    end
    set(gca,'ytick',0:0.05:0.20)
    ylim([-0.01 0.20])
    xlim([0.25 size(temp,2)+0.75])
    box on
    ylabel('Normalized distance')
    axis square
    maxVal = nanmax(temp);
    title(sprintf('%s',runString{iRun}))
    cnt = cnt + 1;
end

% 
subplot(3,3,cnt)
xV = 1;
for iMod = 1:3   
    if iMod == 1
        xV = 1;
    elseif iMod == 2
        xV = 5;
    else
        xV = 9;
    end
    
    for iRun = 6:8
        data = pkDiffNormFinal{iRun,1}(:,iMod);
        data = data(~isnan(data));
        nBoot = 1000;
        
        bootMedians = zeros(nBoot,1);
        for i = 1:nBoot
            resampleData = datasample(data,length(data));
            bootMedians(i) = median(resampleData);
        end
        ci = prctile(bootMedians,[2.5 97.5]);
        med = nanmedian(data);
        
        plot([xV xV],[ci(1) ci(2)],'-','linewidth',1,'color',[0 0.4470 0.7410])
        hold on
        plot(xV,med,'.','markersize',20,'color',[0 0.4470 0.7410])
        
        xV = xV + 1;
    end  
end
set(gca,'xtick',2:4:12,'xticklabels',{'M1','M2','M3'})
set(gca,'ytick',0:0.05:0.20)
xlim([0 12])
ylim([-0.01 0.20])
box on
ylabel('Change in location')
axis square
title('d3')
cnt = cnt + 1;

% 
subplot(3,3,cnt)
xV = 1;
for iMod = 1:3
   
    if iMod == 1
        xV = 1;
    elseif iMod == 2
        xV = 5;
    else
        xV = 9;
    end
    
    for iRun = 6:8
        data = pkDiffNormFinal{iRun,1}(:,iMod);
        data = data(~isnan(data));
        nBoot = 1000;
        
        bootMedians = zeros(nBoot,1);
        for i = 1:nBoot
            resampleData = datasample(data,length(data));
            bootMedians(i) = median(resampleData);
        end
        ci = prctile(bootMedians,[2.5 97.5]);
        med = nanmedian(data);
        
        plot([xV xV],[ci(1) ci(2)],'-','linewidth',1,'color',[0 0.4470 0.7410])
        hold on
        plot(xV,med,'.','markersize',20,'color',[0 0.4470 0.7410])
        
        xV = xV + 1;
    end  
end

set(gca,'xtick',2:4:12,'xticklabels',{'M1','M2','M3'})
set(gca,'ytick',0:0.05:0.20)
xlim([0 12])
ylim([-0.01 0.20])
box on
ylabel('Change in location')
axis square
title('d4')
cnt = cnt + 1;

% d5
subplot(3,3,cnt)
xV = 1;
for iMod = 1:4
   
    if iMod == 1
        xV = 1;
    elseif iMod == 2
        xV = 5;
    elseif iMod == 3
        xV = 9;
    elseif iMod == 4
        xV = 13;
    end
    
    for iRun = 9:11
        data = pkDiffNormFinal{iRun,1}(:,iMod);
        data = data(~isnan(data));
        nBoot = 1000;
        
        bootMedians = zeros(nBoot,1);
        for i = 1:nBoot
            resampleData = datasample(data,length(data));
            bootMedians(i) = median(resampleData);
        end
        ci = prctile(bootMedians,[2.5 97.5]);
        med = nanmedian(data);
        
        plot([xV xV],[ci(1) ci(2)],'-','linewidth',1,'color',[0 0.4470 0.7410])
        hold on
        plot(xV,med,'.','markersize',20,'color',[0 0.4470 0.7410])
        
        xV = xV + 1;
    end  
end

set(gca,'xtick',2:4:16,'xticklabels',{'M1','M2','M2','M3'})
set(gca,'ytick',0:0.05:0.20)
xlim([0 16])
ylim([-0.01 0.20])
box on
ylabel('Change in location')
axis square
title('d5')
cnt = cnt + 1;

% bin size
subplot(3,3,cnt)
xV = 1;
for iMod = 1:3   
    if iMod == 1
        xV = 1;
    elseif iMod == 2
        xV = 6;
    else
        xV = 11;
    end
    
    for iRun = 12:15
        data = pkDiffNormFinal{iRun,1}(:,iMod);
        data = data(~isnan(data));
        nBoot = 1000;
        
        bootMedians = zeros(nBoot,1);
        for i = 1:nBoot
            resampleData = datasample(data,length(data));
            bootMedians(i) = median(resampleData);
        end
        ci = prctile(bootMedians,[2.5 97.5]);
        med = nanmedian(data);
        
        plot([xV xV],[ci(1) ci(2)],'-','linewidth',1,'color',[0 0.4470 0.7410])
        hold on
        plot(xV,med,'.','markersize',20,'color',[0 0.4470 0.7410])
        
        xV = xV + 1;
    end  
end
set(gca,'xtick',2.5:5:15,'xticklabels',{'M1','M2','M3'})
set(gca,'ytick',0:0.05:0.20)
xlim([0 15])
ylim([-0.01 0.20])
box on
ylabel('Change in location')
axis square
title('Bin size')
cnt = cnt + 1;

%% BOOT ROT
runString = {'e1','e2'};

for iRun = 1:2
    if iRun == 1 || iRun == 3
        figure('position',[-1042,15,991,853]);
        cnt = 1;
    end
    
    subplot(3,3,cnt)
    for iMod = 1:3
        data = rotDiffFinal{iRun,1}(:,iMod);
        data = data(~isnan(data));
        nBoot = 1000;
        
        bootMedians = zeros(nBoot,1);
        for i = 1:nBoot
            resampleData = datasample(data,length(data));
            bootMedians(i) = median(resampleData);
        end
        ci = prctile(bootMedians,[2.5 97.5]);
        med = nanmedian(data);

        plot([iMod iMod],[ci(1) ci(2)],'-','linewidth',1,'color',[0 0.4470 0.7410])
        hold on
        plot(iMod,med,'.','markersize',20,'color',[0 0.4470 0.7410])
    end
    
    temp = rotDiffFinal{iRun,1}(:,1:3);
    if size(temp,2) > 3
        set(gca,'xtick',1:size(temp,2),'xticklabels',{'M1','M2','M3','M4'})  
    else
        set(gca,'xtick',1:size(temp,2),'xticklabels',{'M1','M2','M3'})
    end
    set(gca,'ytick',0:15:100)
    ylim([-1 45])
    xlim([0.25 size(temp,2)+0.75])
    box on
    ylabel('Rotation difference')
    axis square
    maxVal = nanmax(temp);
    title(sprintf('%s',runString{iRun}))
    cnt = cnt + 1;
end

% e3
subplot(3,3,cnt)
xV = 1;
for iMod = 1:3   
    if iMod == 1
        xV = 1;
    elseif iMod == 2
        xV = 5;
    else
        xV = 9;
    end
    
    for iRun = 6:8
        data = rotDiffFinal{iRun,1}(:,iMod);
        data = data(~isnan(data));
        nBoot = 1000;
        
        bootMedians = zeros(nBoot,1);
        for i = 1:nBoot
            resampleData = datasample(data,length(data));
            bootMedians(i) = median(resampleData);
        end
        ci = prctile(bootMedians,[2.5 97.5]);
        med = nanmedian(data);
        
        plot([xV xV],[ci(1) ci(2)],'-','linewidth',1,'color',[0 0.4470 0.7410])
        hold on
        plot(xV,med,'.','markersize',20,'color',[0 0.4470 0.7410])
        
        xV = xV + 1;
    end  
end
set(gca,'xtick',2:4:12,'xticklabels',{'M1','M2','M3'})
set(gca,'ytick',0:15:100)
xlim([0 12])
ylim([-1 45])
box on
ylabel('Rotation difference')
axis square
title('e3')
cnt = cnt + 1;

% e4
subplot(3,3,cnt)
xV = 1;
for iMod = 1:3
   
    if iMod == 1
        xV = 1;
    elseif iMod == 2
        xV = 5;
    else
        xV = 9;
    end
    
    for iRun = 6:8
        data = rotDiffFinal{iRun,1}(:,iMod);
        data = data(~isnan(data));
       
        nBoot = 1000;
        
        bootMedians = zeros(nBoot,1);
        for i = 1:nBoot
            resampleData = datasample(data,length(data));
            bootMedians(i) = median(resampleData);
        end
        ci = prctile(bootMedians,[2.5 97.5]);
        med = nanmedian(data);
        
        plot([xV xV],[ci(1) ci(2)],'-','linewidth',1,'color',[0 0.4470 0.7410])
        hold on
        plot(xV,med,'.','markersize',20,'color',[0 0.4470 0.7410])
        
        xV = xV + 1;
    end  
end

set(gca,'xtick',2:4:12,'xticklabels',{'M1','M2','M3'})
set(gca,'ytick',0:15:45)
xlim([0 12])
ylim([-1 45])
box on
ylabel('Rotation difference')
axis square
title('e4')
cnt = cnt + 1;

% e5
subplot(3,3,cnt)
xV = 1;
for iMod = 1:4
   
    if iMod == 1
        xV = 1;
    elseif iMod == 2
        xV = 5;
    elseif iMod == 3
        xV = 9;
    else
        xV = 13;
    end
    
    for iRun = 9:11
        data = rotDiffFinal{iRun,1}(:,iMod);
        data = data(~isnan(data));
       
        nBoot = 1000;
        
        bootMedians = zeros(nBoot,1);
        for i = 1:nBoot
            resampleData = datasample(data,length(data));
            bootMedians(i) = median(resampleData);
        end
        ci = prctile(bootMedians,[2.5 97.5]);
        med = nanmedian(data);
        
        plot([xV xV],[ci(1) ci(2)],'-','linewidth',1,'color',[0 0.4470 0.7410])
        hold on
        plot(xV,med,'.','markersize',20,'color',[0 0.4470 0.7410])
        
        xV = xV + 1;
    end  
end

set(gca,'xtick',2:4:16,'xticklabels',{'M1','M2','M2','M3'})
set(gca,'ytick',0:15:45)
xlim([0 16])
ylim([-1 45])
box on
ylabel('Rotation difference')
axis square
title('e5')
cnt = cnt + 1;

% bin size
subplot(3,3,cnt)
xV = 1;
for iMod = 1:3   
    if iMod == 1
        xV = 1;
    elseif iMod == 2
        xV = 6;
    else
        xV = 11;
    end
    
    for iRun = 12:15
        data = rotDiffFinal{iRun,1}(:,iMod);
        data = data(~isnan(data));
        nBoot = 1000;
        
        bootMedians = zeros(nBoot,1);
        for i = 1:nBoot
            resampleData = datasample(data,length(data));
            bootMedians(i) = median(resampleData);
        end
        ci = prctile(bootMedians,[2.5 97.5]);
        med = nanmedian(data);
        
        plot([xV xV],[ci(1) ci(2)],'-','linewidth',1,'color',[0 0.4470 0.7410])
        hold on
        plot(xV,med,'.','markersize',20,'color',[0 0.4470 0.7410])
        
        xV = xV + 1;
    end  
end
set(gca,'xtick',2.5:5:15,'xticklabels',{'M1','M2','M3'})
set(gca,'ytick',0:15:45)
xlim([0 15])
ylim([-1 45])
box on
ylabel('Rotation difference')
axis square
title('Bin size')
cnt = cnt + 1;

%
subplot(3,3,cnt)
xV = 1;
for iMod = 1:3   
    if iMod == 1
        xV = 1;
    elseif iMod == 2
        xV = 4;
    else
        xV = 7;
    end
    
    for iRun = 16:17
        data = rotDiffFinal{iRun,1}(:,iMod);
        data = data(~isnan(data));
        nBoot = 1000;
        
        bootMedians = zeros(nBoot,1);
        for i = 1:nBoot
            resampleData = datasample(data,length(data));
            bootMedians(i) = median(resampleData);
        end
        ci = prctile(bootMedians,[2.5 97.5]);
        med = nanmedian(data);
        
        plot([xV xV],[ci(1) ci(2)],'-','linewidth',1,'color',[0 0.4470 0.7410])
        hold on
        plot(xV,med,'.','markersize',20,'color',[0 0.4470 0.7410])
        
        xV = xV + 1;
    end  
end
set(gca,'xtick',1.5:3:15,'xticklabels',{'M1','M2','M3'})
set(gca,'ytick',0:15:45)
xlim([0 9])
ylim([-1 45])
box on
ylabel('Rotation difference')
axis square
title('Deg size')
cnt = cnt + 1;
