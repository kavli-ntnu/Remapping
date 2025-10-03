load expStore_ed8c

figure;
cnt = 1;

for iExp = 1:size(expStore,1)
    if ~isempty(expStore{iExp,1})
        
        subplot(2,size(expStore,1)+1,cnt);
        
        pairDistA = expStore{iExp,3};
        pairDistB = expStore{iExp,4};
        
        x = pairDistA;
        y = pairDistB;
        
        % Bin the data:
        % pts = linspace(-4, 4, 101);
        % N = histcounts2(y(:), x(:), pts, pts);
        
        pts = linspace(0, 40*sqrt(2), 80);
        N = histcounts2(y(:), x(:), pts, pts);
        
        % Plot scattered data (for comparison):
        %     figure;
        %     subplot(1, 2, 1);
        %     scatter(x, y, 'r.');
        %     axis equal;
        %     set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]));
        
        % Plot heatmap:
        %     subplot(1, 2, 2);
        imagesc(pts, pts, N);
        axis equal;
        set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
        title(sprintf('%.2f',expStore{iExp,5}))
        cnt = cnt + 1;
    end
end

load expStoreStab_ed8c

for iExp = 1:size(expStoreStab,1)
    if ~isempty(expStoreStab{iExp,1})
        
        subplot(2,size(expStore,1)+1,cnt);
        
        pairDistA = expStoreStab{iExp,3};
        pairDistB = expStoreStab{iExp,4};
        
        x = pairDistA;
        y = pairDistB;
        
        % Bin the data:
        % pts = linspace(-4, 4, 101);
        % N = histcounts2(y(:), x(:), pts, pts);
        
        pts = linspace(0, 40*sqrt(2), 80);
        N = histcounts2(y(:), x(:), pts, pts);
        
        % Plot scattered data (for comparison):
        %     figure;
        %     subplot(1, 2, 1);
        %     scatter(x, y, 'r.');
        %     axis equal;
        %     set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]));
        
        % Plot heatmap:
        %     subplot(1, 2, 2);
        imagesc(pts, pts, N);
        axis equal;
        set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
        title(sprintf('%.2f',expStoreStab{iExp,5}))
        cnt = cnt + 1;
    end
end