load sim5617_hex

figure;
plot(closestM1(:,1),closestM1(:,2),'o')
xline(0)
yline(0)
axis square

hold on
plot(closestM2(:,1),closestM2(:,2),'o')
xline(0)
yline(0)
axis square
  
hold on
plot(closestM3(:,1),closestM3(:,2),'o')
xline(0)
yline(0)
axis square

for iMod = 1:3
    plot(gShift(iMod,1),gShift(iMod,2),'k*')
end

xlim([-110 110])
ylim([-110 110])
title('Sim 5617')

%%
load sim7438_hex

figure;
plot(closestM1(:,1),closestM1(:,2),'o')
xline(0)
yline(0)
axis square

hold on
plot(closestM2(:,1),closestM2(:,2),'o')
xline(0)
yline(0)
axis square
  
hold on
plot(closestM3(:,1),closestM3(:,2),'o')
xline(0)
yline(0)
axis square

for iMod = 1:3
    plot(gShift(iMod,1),gShift(iMod,2),'k*')
end

xlim([-110 110])
ylim([-110 110])
title('Sim 7438')
