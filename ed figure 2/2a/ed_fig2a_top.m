expNums = [1 2];

pltCnt = 1;
fig = figure('position',[1919,-57,1920,963]);

for iExp = expNums
    
    if iExp == 1
        load wkSpace_26821
    elseif iExp == 2
        load wkSpace_27150
        close(gcf)
    end
    
    %% Use DBSCAN to cluster the points in the 2-D UMAP output
    % db_eps = 0.3;     % low values requires dense clusters, high values allow diffuse clusters
    % db_minPts = 10;     % low values allow diffuse clusters, high values require dense clusters
    mcluIds = dbscan(XUmap, db_eps, db_minPts);
    
    % Check how many clusters we have
    mcluIdsU = unique(mcluIds);
    nclu = numel(mcluIdsU);
    fprintf("Found %u clusters\n", nclu);
    
    % % Save results of this umap run for comparison
    % fn = sprintf("umap_res_%s.mat", datestr(now(), "yyyy-mm-dd_HH-MM-SS"));
    % save(fn, "X", "XUmap", "mcluIds", "nclu", "-v7.3", "-nocompression");
    
    % Gather information about each cluster
    clear moduleClu
    
    for c = 1:nclu
        
        clear mclu
        
        % Get indices of all units in the current module cluster
        id = mcluIdsU(c);
        v = mcluIds==id;
        iu = find(v);
        
        % Store unit identity info in module struct
        mclu.id = id;
        mclu.unit_indices = iu;
        tmpIDs = 1:size(XUmap,1);
        mclu.unit_ids = [tmpIDs(v)]';
        mclu.n_units = sum(v);
        %     mclu.acorrs = reshape(acAll(:, v), acsz(1), acsz(2), []);
        
        for u = 1:mclu.n_units
            mclu.ratemaps(:, :, u) = dat.maps{iu(u)};
        end
        
        % Infer whether what cluster type this is (possibilities are "noise",
        % "non-grid", or "grid"). The DBSCAN output is predictable enough to
        % infer this from the cluster ID.
        %
        % "noise":    Always assigned to ID -1, but not present in every
        %             clustering result.
        %
        % "non-grid": In a good outcome, this should be the single largest
        %             cluster. In every such case which I've seen, this cluster
        %             has been assigned to ID 1.
        %
        % "grid":    In a good outcome, all cluster IDs greater than 1 should
        %            correspond to grid modules.
        
        % Ascertain the type of cluster, based on the rules/assumptions above
        if mclu.id == -1
            mclu.name = "noise";
        elseif mclu.id == 1
            mclu.name = "non-grid";
        elseif mclu.id > 1
            mclu.name = "grid";
        end
        moduleClu(c) = mclu;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % To make things easier, re-order the grid clusters according to their
    % % grid spacing
    %
    gscores = dat.grids(:,1);
    meanGspacing = dat.grids(:,2);
    isgridclu = [moduleClu.id] >= 2;
    moduleCluGrid = moduleClu(isgridclu);
    ngridclu = sum(isgridclu);
    clear cluGridSpacing
    cluGridSpacing = [];
    for c = 1:ngridclu
        mclu = moduleCluGrid(c);
        inds = mclu.unit_indices;
        cluGridSpacing(c) = nanmean(meanGspacing(inds));
    end
    
    % sort grid clusters by spacing and rename
    [cluGridSpacing, isort] = sort(cluGridSpacing,'descend');
    moduleCluGrid = moduleCluGrid(isort);
    for c = 1:ngridclu
        moduleCluGrid(c).name = sprintf("grid M%u", c);
        moduleCluGrid(c).id = c+1;
    end
    moduleClu = [moduleClu(~isgridclu), moduleCluGrid];
    
    mcluIds = nan(nu, 1);
    for c = 1:nclu
        mclu = moduleClu(c);
        mcluIds(mclu.unit_indices) = mclu.id;
        fprintf("cluster ID %d, '%s', n=%u units\n", mclu.id, mclu.name, mclu.n_units);
    end
    disp(' ');
    
    % PLOT CLUSTERING RESULTS
    
    if size(expNums,2) == 4
        subplot(2,2,pltCnt)
    else
        subplot(3,1,pltCnt)
    end
    
    hold on
    clear cols
    
    cols(1, :) = 0.5 + [0, 0, 0];
    cols(2, :) = [0.8, 0, 0];
    cols(3, :) = [0.8, 0.6, 0.0];
    cols(4, :) = [0, 0.3, 1.0];
    cols(5, :) = [0, 0.7, 0.2];
    
    cols = [cols; colorcube(20)];
    
    clear h igrid ngrid
    cnt = 0;
    x = XUmap(:, 1);
    y = XUmap(:, 2);
    if size(XUmap, 2) > 2
        z = XUmap(:, 3);
    else
        z = zeros(nu, 1);
    end
    
    for c = 1:nclu
        mclu = moduleClu(c);
        
        vclu = mcluIds==mclu.id;
        if mclu.name == "noise"
            h(c) = plot(x(vclu), y(vclu), "kx");
        else
            cnt = cnt+1;
            h(c) = plot(x(vclu), y(vclu), ".", "color", cols(cnt, :));
        end
        %     try
        %         dtip = dataTipTextRow("id", mclu.unit_ids);
        %         h(c).DataTipTemplate.DataTipRows(end+1) = dtip;
        %     end
        legstrs(c) = sprintf("%-10sn = %u", mclu.name, sum(vclu));
    end
    
    colormap lines
    % [~,icons] = legend(h, legstrs, "fontName", "consolas", "position", [0.7264    0.8448    0.2515    0.1170]);
%     [~,icons] = legend(h, legstrs, "fontName", "consolas", "location", "northeastoutside");
    
    set(findobj(icons,'type','line'),'markers',25)
    xlabel 'UMAP dim 1'
    ylabel 'UMAP dim 2'
    axis square
    
    pltCnt = pltCnt + 1;
end
