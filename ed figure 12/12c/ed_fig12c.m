%%
load simDist_simDisp_ed12c

figure;
[f1,x1] = ksdensity(disp(:,1),0:0.01:tan(pi/6),'bandwidth',0.1);
[pks,loc,width,prom] = findpeaks(f1,x1)
findpeaks(f1,x1)
hold on

[f2,x2] = ksdensity(disp(:,2),0:0.01:tan(pi/6),'bandwidth',0.1);
[pks,loc,width,prom] = findpeaks(f2,x2)
findpeaks(f2,x2)
hold on

[f3,x3] = ksdensity(disp(:,3),0:0.01:tan(pi/6),'bandwidth',0.1);
[pks,loc,width,prom] = findpeaks(f3,x3)
findpeaks(f3,x3)
hold on

xlabel('Normalized displacement')
ylabel('Proportion')

%%
figure;
[f1,x1] = ksdensity(dist(:,1),0:0.01:1,'bandwidth',0.1);
[pks,loc,width,prom] = findpeaks(f1,x1)
findpeaks(f1,x1)
hold on

[f2,x2] = ksdensity(dist(:,2),0:0.01:1,'bandwidth',0.1);
[pks,loc,width,prom] = findpeaks(f2,x2)
findpeaks(f2,x2)
hold on

[f3,x3] = ksdensity(dist(:,3),0:0.01:1,'bandwidth',0.1);
[pks,loc,width,prom] = findpeaks(f3,x3)
findpeaks(f3,x3)
hold on

ylim([0 2.2])
xlabel('Normalized distance')
ylabel('Proportion')