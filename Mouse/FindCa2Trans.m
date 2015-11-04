function [ output_args ] = FindCa2Trans( transData )

[~,locs]=findpeaks(smooth(smooth(transData)),'MINPEAKHEIGHT', 30);

for i=1:length(locs)
    
    newData(i,:) = transData([locs(i)-100:locs(i)+100]);
    
end % END FOR

meanLine = mean(newData,1);

errorData = std(newData,0,1);

meanLine = smooth(meanLine);
errorData = smooth(errorData);

patchTime = [time, fliplr(time)];

patchData = smooth([meanLine'+errorData', fliplr(meanLine'-errorData')]);

figure(2)
hold on
lData = plot(time(1:2:end), meanLine(1:2:end), 'r', 'linewidth', 1);
pData = patch(patchTime, patchData, 1);
lData = plot(time(1:5:end), meanLine(1:5:end), 'r', 'linewidth', 2.75);
hold off

set(pData, 'FaceColor', 'r')
set(pData, 'EdgeColor', 'none')
set(pData, 'FaceAlpha', 0.25)


end % END FUNCTION

% EOF