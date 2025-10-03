load expStore_29730

figure;
subplot(122)

pairDistA = expStore{1,3};
pairDistB = expStore{1,4};

x = pairDistA;
y = pairDistB;

pts = linspace(0, 40*sqrt(2), 80);
N = histcounts2(y(:), x(:), pts, pts);

imagesc(pts, pts, N);
axis equal;
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
title(sprintf('%.2f',expStore{1,5}))

%%
load expStoreStab_29730
subplot(121)

pairDistA = expStoreStab{1,3};
pairDistB = expStoreStab{1,4};

x = pairDistA;
y = pairDistB;

pts = linspace(0, 40*sqrt(2), 80);
N = histcounts2(y(:), x(:), pts, pts);

imagesc(pts, pts, N);
axis equal;
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
title(sprintf('%.2f',expStoreStab{1,5}))

