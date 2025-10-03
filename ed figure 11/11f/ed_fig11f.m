load hdCells_ed11f

% plot tuning curves of HD cells
for iCell = 1:size(hdCells,1)
    figure;
    hold on
    for iSession = 1:4
        hd = sData{1,iSession}.results.hd;
        spkHD = hd(sData{1,iSession}.results.spkInd{hdCells(iCell)});
        if ~isempty(spkHD)
            tc = analyses.turningCurve(spkHD,hd,0.0083,'binWidth',6);
            tcStat = analyses.tcStatistics(tc,6,20);
            
            subplot(1,4,iSession)
            circularTurningBRK(tc(:,2)/max(tc(:,2)),'k-','linewidth',3)
            circularTurningBRK(tc(:,3)/max(tc(:,3)),'adjustaxis',false,'color',[.5 .5 .5])
            axis equal
            title(sprintf('length = %.2f angle = %.2f',tcStat.r,mod(360-tcStat.mean,360)),'fontweight','normal','fontsize',10);
        end
    end
    
end