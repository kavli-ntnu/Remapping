load fxCells_27207

load colorMap
xVals = [1 1.5];
selectedCol = 3;

figure('position',[1295,238,1218,374])
subplot(141)

numMods = 2;
for iMod = 1:numMods
    gridDepth = gridCells(gridCells(:,2) == iMod,5);
    plotSpread(gridDepth,'xValues',xVals(iMod))
    o = flipud(findobj(gca,'type','line'));
    set(o(iMod),'color',colorMap(iMod,:))
    set(gca,'fontsize',12)
end

hdDepth = nan(size(hdCells,1),1);
for iRow = 1:size(hdCells,1)   
    hdDepth(iRow,1) = 1;
    hdDepth(iRow,2) = hdCells(iRow,end);
    hdDepth(iRow,3) = diffMeanAngle(iRow,selectedCol);
end
hd1 = hdDepth(hdDepth(:,1) == 1,:);

plotSpread(hd1(:,2),'xValues',numMods+1)
p = flipud(findobj(gca,'type','line'));
set(p(size(p,1)),'color',colorMap(6,:))

if ~isempty(borderCells)
    plotSpread(borderCells(:,end),'xValues',numMods+2)
    p = flipud(findobj(gca,'type','line'));
    set(p(size(p,1)),'color',colorMap(9,:))
end

ylabel('Distance from tip (um)')
hold on
if numMods > 0
    set(gca,'xtick',[nanmedian(xVals) numMods+1 numMods+2],'xticklabels',{'Grid','HD','Border'})
else
    set(gca,'xtick',[1 2],'xticklabels',{'HD','Border'})
end

xtickangle(45)

%
subplot(142)
plotSpread(gridCellRot)
p = flipud(findobj(gca,'type','line'));

hold on

for iMod = 1:numMods
    plot(iMod,nanmedian(gridCellRot{1,iMod}),'r+','markersize',15,'linewidth',2)
end
ylabel('Rotation (deg)')
xlim([0.5 2.5])
ylim([0 360])
box off
set(gca,'xtick',1:2,'xticklabels',{'M1','M2'})
set(gca,'fontsize',12)

%
subplot(143)
plotSpread(hd1(:,3),'xValues',0.8)
p = flipud(findobj(gca,'type','line'));
set(p(1),'color',colorMap(6,:))
hold on

HDmediansXdepth(1) = meanangle(hd1(:,3),[0 360]);
plot(0.8,HDmediansXdepth(1),'x','color','k','markersize',10,'linewidth',2)

set(gca,'xtick','')
xlabel('HD')
ylim([0 360])
set(gca,'fontsize',12)


subplot(144)
plotSpread(borderRotAngle)
p = flipud(findobj(gca,'type','line'));
set(p(1),'color',colorMap(9,:))

box off
xlabel('Border')
ylim([0 360])
set(gca,'fontsize',12)
set(gca,'xtick','')

