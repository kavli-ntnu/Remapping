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

%
spacing = [50 70 100];
halfSpacing = spacing * 0.5;

colorblindMap
% figure;
for iMod = 1:3
    hexCoords = [halfSpacing(iMod) 0; halfSpacing(iMod) halfSpacing(iMod)*tan(pi/6); 0 halfSpacing(iMod); -halfSpacing(iMod) halfSpacing(iMod)*tan(pi/6); -halfSpacing(iMod) 0; -halfSpacing(iMod) -halfSpacing(iMod)*tan(pi/6); 0 -halfSpacing(iMod); halfSpacing(iMod) -halfSpacing(iMod)*tan(pi/6); halfSpacing(iMod) 0];
    plot(hexCoords(:,1),hexCoords(:,2),'k-')
    hold on
    %     hexStore{1,iMod} = hexCoords;
end
xline(0)
yline(0)
axis square
xlim([-205 205])
ylim([-205 205])
% set(gca,'xtick',-50:25:50,'ytick',-50:25:50)

boxCoords = [-100 -100; 100 -100; 100 100; -100 100; -100 -100];
plot(boxCoords(:,1),boxCoords(:,2),'k-')
axis square
close(gcf)

%%
selectedSims = [3202 4283 390];
load rotationsToTest

PLOT = true;
numRuns = 1;

for iSim = 1:size(rotationsToTest,1)
    
    % apply rotation
    rotAngle = rotationsToTest{iSim,1};
    
    rotM1(:,1) =  m1(:,1)*cosd(rotAngle(1)) + m1(:,2)*sind(rotAngle(1));
    rotM1(:,2) = -m1(:,1)*sind(rotAngle(1)) + m1(:,2)*cosd(rotAngle(1));
    
    rotM2(:,1) =  m2(:,1)*cosd(rotAngle(2)) + m2(:,2)*sind(rotAngle(2));
    rotM2(:,2) = -m2(:,1)*sind(rotAngle(2)) + m2(:,2)*cosd(rotAngle(2));
    
    rotM3(:,1) =  m3(:,1)*cosd(rotAngle(3)) + m3(:,2)*sind(rotAngle(3));
    rotM3(:,2) = -m3(:,1)*sind(rotAngle(3)) + m3(:,2)*cosd(rotAngle(3));
    
    if PLOT == true
        figure;
        plot(rotM1(:,1),rotM1(:,2),'ko')
        hold on
        plot(rotM2(:,1),rotM2(:,2),'ro')
        plot(rotM3(:,1),rotM3(:,2),'bo')
        
        boxCoords = [-100 -100; 100 -100; 100 100; -100 100; -100 -100];
        plot(boxCoords(:,1),boxCoords(:,2),'k-')
        axis square
        xline(0); yline(0)
        xlim([-105 105])
        ylim([-105 105])
        title(sprintf('Sim %d',selectedSims(1,iSim)))
    end
end
