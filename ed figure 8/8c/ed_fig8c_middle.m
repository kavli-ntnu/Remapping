%%
load pkStore_ed8c

figure;
pltCnt = 1;

for iRow = 1:size(expStore,1)
    numCells = size(expStore{iRow,3},1);
    colorPlt = 1; % 1 for colored dots, 0 for black arrows only
    
    a = rand([numCells 1]);
    b = rand([numCells 1]);
    c = rand([numCells 1]);
    
    colors = horzcat(a,b,c);
    cnt = 1;
    
    pkStore = expStore{iRow,3};
    
    for iCell = 1:numCells
        c1 = pkStore(iCell,1);
        r1 = pkStore(iCell,2);
        c2 = pkStore(iCell,3);
        r2 = pkStore(iCell,4);
        
        if ~isempty(c1) & ~isempty(r1) & ~isempty(c2) & ~isempty(r2)
            
            subplot(4,2,pltCnt)
            if colorPlt == 1
                plot(c1(1),r1(1),'.','markersize',35,'color',colors(cnt,:))
                hold on
                plot(c2(1),r2(1),'.','markersize',35,'color',colors(cnt,:))
            end
            arrow([c1(1),r1(1)],[c2(1) r2(1)],11)
            axis square
            xlim([0 40])
            ylim([0 40])
            
            p1 = [c1(1) r1(1)]; % First Point
            p2 = [c2(1) r2(1)];
            dp = p2-p1; % Difference
            % quiver(p1(1),p1(2),dp(1),dp(2),0,'markersize',50)
            
            cnt = cnt + 1;
            set(gca,'xticklabels','','yticklabels','')
            set(gca,'ydir','rev')
            box on
        end
    end
    pltCnt = pltCnt + 1;
end


%%
load pkStoreStab_ed8c
expStore = expStoreStab;

figure;
pltCnt = 1;

for iRow = 1 %:size(expStore,1)
    numCells = size(expStore{iRow,3},1);
    colorPlt = 1; % 1 for colored dots, 0 for black arrows only
    
    a = rand([numCells 1]);
    b = rand([numCells 1]);
    c = rand([numCells 1]);
    
    colors = horzcat(a,b,c);
    cnt = 1;
    
    pkStore = expStore{iRow,3};
    
    for iCell = 1:numCells
        c1 = pkStore(iCell,1);
        r1 = pkStore(iCell,2);
        c2 = pkStore(iCell,3);
        r2 = pkStore(iCell,4);
        
        if ~isempty(c1) & ~isempty(r1) & ~isempty(c2) & ~isempty(r2)
            
            subplot(4,2,pltCnt)
            if colorPlt == 1
                plot(c1(1),r1(1),'.','markersize',35,'color',colors(cnt,:))
                hold on
                plot(c2(1),r2(1),'.','markersize',35,'color',colors(cnt,:))
            end
            arrow([c1(1),r1(1)],[c2(1) r2(1)],11)
            axis square
            xlim([0 40])
            ylim([0 40])
            
            p1 = [c1(1) r1(1)]; % First Point
            p2 = [c2(1) r2(1)];
            dp = p2-p1; % Difference
            % quiver(p1(1),p1(2),dp(1),dp(2),0,'markersize',50)
            
            cnt = cnt + 1;
            set(gca,'xticklabels','','yticklabels','')
            set(gca,'ydir','rev')
            box on
        end
    end
    pltCnt = pltCnt + 1;
end


