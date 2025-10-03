load CA3_sessionA

mapsBL = p3R;
validBL = i3R;

%%
simsToRun = [1382 4858 7347];

for iSim = 1:size(simsToRun,2)
    simNum = simsToRun(iSim);
    
    if iSim == 1
        load CA3_sessionB_1382
    elseif iSim == 2
        load CA3_sessionB_4858
    else
        load CA3_sessionB_7347
    end

    maps1 = p3R;
    valid1 = i3R;
    
    %
    figure('position',[1,41,1280,683]);
    numCells = sum(validBL == 1 & valid1 == 1)
    
    a = rand([numCells 1]);
    b = rand([numCells 1]);
    c = rand([numCells 1]);
    
    colors = horzcat(a,b,c);
    
    scTest = nan(size(validBL,1),1);
    dpStore = nan(size(validBL,1),2);
    cnt = 1;
    for iCell = 1:size(validBL,1)
        if validBL(iCell) && valid1(iCell)
            
            subplot(131)
            mapBL = general.smooth(mapsBL(:,:,iCell),6.67);
            
            [r1,c1] = find(mapBL == nanmax(nanmax(mapBL)));
            
            plot(c1(1),r1(1),'.','markersize',35,'color',colors(cnt,:))
            hold on
            axis square
            
            subplot(132)
            
            mapCNO = general.smooth(maps1(:,:,iCell),6.67);
            %         [r2,c2] = find(mapCNO == nanmax(nanmax(mapCNO)));
            
            [fmap,fields] = analyses.placefield(mapCNO);
            if ~isempty(fields)
                fieldLoc = nan(size(fields,2),2);
                for iField = 1:size(fields,2)
                    fieldLoc(iField,1) = fields(iField).peakX; %c
                    fieldLoc(iField,2) = fields(iField).peakY; %r
                end
                [idx,d] = knnsearch(fieldLoc,[c1(1) r1(1)]);
                
                dStore(iCell,1) = d;
                
                c2 = fieldLoc(idx,1);
                r2 = fieldLoc(idx,2);
                
                plot(c2(1),r2(1),'.','markersize',35,'color',colors(cnt,:))
                hold on
                axis square
                
                subplot(133)
                plot(c1(1),r1(1),'.','markersize',35,'color',colors(cnt,:))
                hold on
                plot(c2(1),r2(1),'.','markersize',35,'color',colors(cnt,:))
                arrow([c1(1),r1(1)],[c2(1) r2(1)],5)
                axis square
                xlim([0 100])
                ylim([0 100])
                
                p1 = [c1(1) r1(1)]; % First Point
                p2 = [c2(1) r2(1)];
                dp = p2-p1; % Difference
                %     quiver(p1(1),p1(2),dp(1),dp(2),0,'markersize',50)
                
                dpStore(iCell,:) = dp;
                
                scTest(iCell,1) = analyses.spatialCrossCorrelation(mapBL,mapCNO);
                cnt = cnt + 1;
                set(gca,'xticklabels','','yticklabels','')
            end
        end
    end
end
