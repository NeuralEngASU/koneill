% Plot the calcium transients from GenCa2Trans


function [  ] = PlotCa2Trans( transData, time )

meanLine = mean(transData,1);

errorData = std(transData,0,1);

meanLine = smooth(meanLine);
errorData = smooth(errorData);

patchTime = [time, fliplr(time)];

patchData = smooth([meanLine'+errorData', fliplr(meanLine'-errorData')]);

figure(2)
hold on
lData = plot(time(1:5:end), meanLine(1:5:end), 'r', 'linewidth', 1);
pData = patch(patchTime, patchData, 1);
lData = plot(time(1:5:end), meanLine(1:5:end), 'r', 'linewidth', 2.75);
hold off

set(pData, 'FaceColor', 'r')
set(pData, 'EdgeColor', 'none')
set(pData, 'FaceAlpha', 0.25)

% figure(2)
% tmpData = [transData(1,:), transData(2,:), transData(3,:), transData(4,:), transData(5,:)];
% plot(0:5/100:24.98,smooth(tmpData(1:5:end)), 'k', 'linewidth', 2.75);

title('Fake Ca^{2+} Transient Data')
xlabel('Time, sec')
ylabel('\DeltaF/F')

set(gcf, 'Units', 'inches')
% set(gcf, 'Position',[0 0 5 1.75])
set(gcf, 'Position',[0 0 5 3])
set(gcf, 'PaperUnits','inches','PaperPosition',[2 2 5 3])

end

