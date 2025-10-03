load finalDistStore_29730
randDist = nanmean(finalDistStore{1,1},2)
figure;
[f,xi] = ksdensity(randDist,'bandwidth',6);
findpeaks(f,xi)
[pk,loc] = findpeaks(f,xi);
peak = loc

xlim([-2 72])
ylim([-0.001 0.04])

