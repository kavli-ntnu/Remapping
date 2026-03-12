load expList_MEC_4e
load fxCellRotation_4e

%%
% module roatation
gridRotation = expStore(:,1);

% HD rotation
hdRot = expStore(:,2);

% border cell rotation
borderRot = expStore(:,3);

%%
figure('position',[723,214,274,314]);
for iExp = 1:size(expList_MEC,1) 
    if sum(~isnan(expStore(iExp,:))) > 1
        xVals = find(~isnan(expStore(iExp,:)));
        yVals = expStore(iExp,~isnan(expStore(iExp,:)));
        plot(xVals,yVals,'-','color',[0.8 0.8 0.8],'linewidth',1.5)
        hold on
        plot(xVals,yVals,'o','color',[0.5 0.5 0.5],'linewidth',1.5)
    end
end
xlim([0.5 3.5])
ylim([-10 375])
set(gca,'ytick',0:60:360)
set(gca,'xtick',1:3,'xticklabels',{'Grid','HD','Border'})
ylabel('Rotation (deg)')
box on
