
%% Glut v TTX
if exist('glutTransMean') && exist('ttxTransMean')
    
    figure(30)
    hold on
    
    % plots for legend

    plot(1:7, glutTransMean(1:end-1), '-or', 'MarkerSize', 3, 'MarkerFaceColor', 'r')
    plot(1:7, ttxTransMean, '-og', 'MarkerSize', 3, 'MarkerFaceColor', 'g')

    
    
        % Plot Glut
    plot(1:7, glutTransMean(1:end-1), '-or', 'MarkerSize', 3, 'MarkerFaceColor', 'r')
    for i = 1:7
        plot([i, i], [glutTransMean(i) + glutErrMat(i), glutTransMean(i) - glutErrMat(i)], '-r')
    end % END FOR
    
    % Plot Glut + TTX
    plot(1:7, ttxTransMean, '-og', 'MarkerSize', 3, 'MarkerFaceColor', 'g')
    for i = 1:7
        plot([i, i], [ttxTransMean(i) + ttxErrMat(i), ttxTransMean(i) - ttxErrMat(i)], '-g')
    end % END FOR
    

    

    
    hold off
    
    title(sprintf('Uncaged Intensity Peak Comparison\n(Glut: n = %d) (Glut+TTX: n = %d)', glutNumTrials, ttxNumTrials ))
    ylabel('%\DeltaF/F')
    xlabel('Uncage Event')
    
    legend({'Glut', 'Glut + TTX'}, 'Location', 'Northwest')
    
    ylim([0, 5])
    xlim([1, 7])
    
    set(gca, 'YTick', [-1, 0, 1, 2, 3, 4, 5, 6])
    set(gca, 'YTickLabel', [-100, 0, 100, 200, 300, 400, 500, 600]);
    
    set(gca, 'XTick', [1:8])
    set(gca, 'XTickLabel', [1:8]);
    
    set(gcf, 'Units', 'inches')
    set(gcf, 'Position', [1 1 3 2.5])
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3 2.5])
    
    set(gca,'TickDir','out')
    
    print('-dpng', ['C:\CodeRepo\Lab\Mouse\Mouse\Figures\GlutvsTTX_AvgStdErr_ErrorBar_UncageStep_All.png'], '-r100');
    
end

%% Glut v DHPG

if exist('glutTransMean') && exist('dhpgTransMean')
    
    figure(31)
    hold on
    
    % plots for legend

    plot(1:7, glutTransMean(1:end-1), '-or', 'MarkerSize', 3, 'MarkerFaceColor', 'r')
    plot(1:7, dhpgTransMean, '-og', 'MarkerSize', 3, 'MarkerFaceColor', 'g')

    
    % Plot Glut
    plot(1:7, glutTransMean(1:end-1), '-or', 'MarkerSize', 3, 'MarkerFaceColor', 'r')
    for i = 1:7
        plot([i, i], [glutTransMean(i) + glutErrMat(i), glutTransMean(i) - glutErrMat(i)], '-r')
    end % END FOR
    
    % Plot Glut + TTX
    plot(1:7, dhpgTransMean, '-og', 'MarkerSize', 3, 'MarkerFaceColor', 'g')
    for i = 1:7
        plot([i, i], [dhpgTransMean(i) + dhpgErrMat(i), dhpgTransMean(i) - dhpgErrMat(i)], '-g')
    end % END FOR
    

    

    
    hold off
    
    title(sprintf('Uncaged Intensity Peak Comparison\n(Glut: n = %d) (DHPG: n = %d)', glutNumTrials, dhpgNumTrials ))
    ylabel('%\DeltaF/F')
    xlabel('Uncage Event')
    
    legend({'Glut', 'DHPG'}, 'Location', 'Northwest')
    
    ylim([0, 5])
    xlim([1, 7])
    
    set(gca, 'YTick', [-1, 0, 1, 2, 3, 4, 5, 6])
    set(gca, 'YTickLabel', [-100, 0, 100, 200, 300, 400, 500, 600]);
    
    set(gca, 'XTick', [1:8])
    set(gca, 'XTickLabel', [1:8]);
    
    set(gcf, 'Units', 'inches')
    set(gcf, 'Position', [1 1 3 2.5])
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3 2.5])
    
    set(gca,'TickDir','out')
    
    print('-dpng', ['C:\CodeRepo\Lab\Mouse\Mouse\Figures\GlutvsDHPG_AvgStdErr_ErrorBar_UncageStep_All.png'], '-r100');
    
end

%% Glut v GABA

if exist('glutTransMean') && exist('gabaTransMean')
    
    figure(32)
    hold on
    
    % plots for legend

    plot(1:7, glutTransMean(1:end-1), '-or', 'MarkerSize', 3, 'MarkerFaceColor', 'r')
    plot(1:7, gabaTransMean, '-og', 'MarkerSize', 3, 'MarkerFaceColor', 'g')

    
    
        % Plot Glut
    plot(1:7, glutTransMean(1:end-1), '-or', 'MarkerSize', 3, 'MarkerFaceColor', 'r')
    for i = 1:7
        plot([i, i], [glutTransMean(i) + glutErrMat(i), glutTransMean(i) - glutErrMat(i)], '-r')
    end % END FOR
    
    % Plot Glut + TTX
    plot(1:7, gabaTransMean, '-og', 'MarkerSize', 3, 'MarkerFaceColor', 'g')
    for i = 1:7
        plot([i, i], [gabaTransMean(i) + gabaErrMat(i), gabaTransMean(i) - gabaErrMat(i)], '-g')
    end % END FOR
    

    

    
    hold off
    
    title(sprintf('Uncaged Intensity Peak Comparison\n(Glut: n = %d) (GABA: n = %d)', glutNumTrials, gabaNumTrials ))
    ylabel('%\DeltaF/F')
    xlabel('Uncage Event')
    
    legend({'Glut', 'GABA'}, 'Location', 'Northwest')
    
    ylim([0, 5])
    xlim([1, 7])
    
    set(gca, 'YTick', [-1, 0, 1, 2, 3, 4, 5, 6])
    set(gca, 'YTickLabel', [-100, 0, 100, 200, 300, 400, 500, 600]);
    
    set(gca, 'XTick', [1:8])
    set(gca, 'XTickLabel', [1:8]);
    
    set(gcf, 'Units', 'inches')
    set(gcf, 'Position', [1 1 3 2.5])
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3 2.5])
    
    set(gca,'TickDir','out')
    
    print('-dpng', ['C:\CodeRepo\Lab\Mouse\Mouse\Figures\GlutvsGABA_AvgStdErr_ErrorBar_UncageStep_All.png'], '-r100');
    
end

%% Glut v TTX v GABA

if exist('glutTransMean') && exist('ttxTransMean') && exist('gabaTransMean')
    
    figure(33)
    hold on
    
    % plots for legend

    plot(1:7, glutTransMean(1:end-1), '-or', 'MarkerSize', 3, 'MarkerFaceColor', 'r')
    plot(1:7, ttxTransMean, '-og', 'MarkerSize', 3, 'MarkerFaceColor', 'g')
    plot(1:7, gabaTransMean, '-ob', 'MarkerSize', 3, 'MarkerFaceColor', 'b')

    
    
        % Plot Glut
    plot(1:7, glutTransMean(1:end-1), '-or', 'MarkerSize', 3, 'MarkerFaceColor', 'r')
    for i = 1:7
        plot([i, i], [glutTransMean(i) + glutErrMat(i), glutTransMean(i) - glutErrMat(i)], '-r')
    end % END FOR
    
    % Plot Glut + TTX
    plot(1:7, ttxTransMean, '-og', 'MarkerSize', 3, 'MarkerFaceColor', 'g')
    for i = 1:7
        plot([i, i], [ttxTransMean(i) + ttxErrMat(i), ttxTransMean(i) - ttxErrMat(i)], '-g')
    end % END FOR
    
    % Plot GABA
    plot(1:7, gabaTransMean, '-ob', 'MarkerSize', 3, 'MarkerFaceColor', 'b')
    for i = 1:7
        plot([i, i], [gabaTransMean(i) + gabaErrMat(i), gabaTransMean(i) - gabaErrMat(i)], '-b')
    end % END FOR
    

    

    
    hold off
    
    title(sprintf('Uncaged Intensity Peak Comparison\n(Glut: n = %d) (Glut+TTX: n = %d) (GABA: n = %d)', glutNumTrials, ttxNumTrials, gabaNumTrials ))
    ylabel('%\DeltaF/F')
    xlabel('Uncage Event')
    
    legend({'Glut', 'Glut + TTX', 'GABA'}, 'Location', 'Northwest')
    
    ylim([0, 5])
    xlim([1, 7])
    
    set(gca, 'YTick', [-1, 0, 1, 2, 3, 4, 5, 6])
    set(gca, 'YTickLabel', [-100, 0, 100, 200, 300, 400, 500, 600]);
    
    set(gca, 'XTick', [1:8])
    set(gca, 'XTickLabel', [1:8]);
    
    set(gcf, 'Units', 'inches')
    set(gcf, 'Position', [1 1 3.5 2.5])
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.5 2.5])
    
    set(gca,'TickDir','out')
    
    print('-dpng', ['C:\CodeRepo\Lab\Mouse\Mouse\Figures\GlutvsTTXvsGABA_AvgStdErr_ErrorBar_UncageStep_All.png'], '-r100');
    
end
