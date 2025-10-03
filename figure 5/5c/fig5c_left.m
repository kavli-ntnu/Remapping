
%% generate the patterns
xVals1 = -200:50:200;
coords = [];
figure;
for y = -173.2051:100*sin(deg2rad(60)):173.2051
    yVals = repmat(y,[1 size(xVals1,2)]);
    plot(xVals1,yVals,'ko')
    hold on
    
    joined = [xVals1' yVals'];
    coords = [coords; joined];
end
% xlim([-200 200])
% ylim([-200 200])

xVals2 = -175:50:175;
for y = -129.9038:100*sin(deg2rad(60)):129.9039
    yVals = repmat(y,[1 size(xVals2,2)]);
    plot(xVals2,yVals,'ko')
    hold on
    
    joined = [xVals2' yVals'];
    coords = [coords; joined];
end

m1 = coords;

xVals1 = -210:70:210;
coords = [];
hold on
for y = -242.4871:140*sin(deg2rad(60)):242.4871
    yVals = repmat(y,[1 size(xVals1,2)]);
    plot(xVals1,yVals,'ro')
    hold on
    
    joined = [xVals1' yVals'];
    coords = [coords; joined];
end

xVals2 = -175:70:175;
for y = -181.8653:140*sin(deg2rad(60)):181.8654
    yVals = repmat(y,[1 size(xVals2,2)]);
    plot(xVals2,yVals,'ro')
    hold on
    
    joined = [xVals2' yVals'];
    coords = [coords; joined];
end

m2 = coords;


xVals1 = -200:100:200;
coords = [];
for y = -173.2051:200*sin(deg2rad(60)):173.2051
    yVals = repmat(y,[1 size(xVals1,2)]);
    plot(xVals1,yVals,'bo')
    hold on
    
    joined = [xVals1' yVals'];
    coords = [coords; joined];
end

xVals2 = -150:100:150;
for y = -259.8076:200*sin(deg2rad(60)):259.8076
    yVals = repmat(y,[1 size(xVals2,2)]);
    plot(xVals2,yVals,'bo')
    hold on
    
    joined = [xVals2' yVals'];
    coords = [coords; joined];
end
xlim([-220 220])
ylim([-220 220])

m3 = coords;

% add the box
boxCoords = [-100 -100; 100 -100; 100 100; -100 100; -100 -100];
plot(boxCoords(:,1),boxCoords(:,2),'k-')
axis square

close(gcf)


%%
load shiftLoc

%% M1
[idx,d] = knnsearch(m1,[0 0],'k',7);
ccCoordStore1 = m1(idx(2:end),:);
ccCoordStore1(end+1,:) = ccCoordStore1(1,:);

x = ccCoordStore1(:,1);
y = ccCoordStore1(:,2);
k = convhull(x,y);

figure; 
plot(x(k),y(k),'b-')
hold on

% M2
[idx,d] = knnsearch(m2,[0 0],'k',7);
ccCoordStore2 = m2(idx(2:end),:);
ccCoordStore2(end+1,:) = ccCoordStore2(1,:);

x = ccCoordStore2(:,1);
y = ccCoordStore2(:,2);
k = convhull(x,y);

plot(x(k),y(k),'r-')

% M3
[idx,d] = knnsearch(m3,[0 0],'k',7);
ccCoordStore3 = m3(idx(2:end),:);
ccCoordStore3(end+1,:) = ccCoordStore3(1,:);

x = ccCoordStore3(:,1);
y = ccCoordStore3(:,2);
k = convhull(x,y);

plot(x(k),y(k),'k-')
axis square
set(gca,'xtick','','ytick','')
xlim([-105 105])
ylim([-105 105])
xline(0)
yline(0)

close(gcf)

%% shifted points
load colorMap
figure;
for i = 1:3
    s = shiftLoc{i,1};
    
    subplot(1,3,i)
    plot([0 s(1,1)],[0 s(1,2)],'-','color',colorMap(1,:),'linewidth',2)
    hold on
    plot([0 s(2,1)],[0 s(2,2)],'-','color',colorMap(2,:),'linewidth',2)
    plot([0 s(3,1)],[0 s(3,2)],'-','color',colorMap(3,:),'linewidth',2)
    
    axis square
    xline(0)
    yline(0)
    xlim([-80 80])
    ylim([-80 80])
end
