load silCoef_ed2b

binRange = -1:0.05:1;
figure; 
histogram(nonGridCol,binRange)
hold on
histogram(finalGridCol,binRange)
ylim([0 42])
xlim([-1 1.05])
set(gca,'ytick',0:10:40)
ylabel('Number of clusters')
xlabel('Silhouette coefficient')
