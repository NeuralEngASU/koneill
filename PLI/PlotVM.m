%% Plot deltaPhi (use while debugging)
figure;
for ii = 1:size(deltaPhi,1)%1907%1047%1700:2100%6050:6150%size(tmpDeltaPhi, 2)
    
%     subplot(2,1,1)
    binEdge = [-pi:pi/100:pi];
    [counts, centers] = hist(squeeze(deltaPhi(ii,1,:)), binEdge);
    P = polar(centers, 122 * ones(1,size(centers,2)));
%     set(P, 'Visible', 'off')
    hold on
    
    polar(centers, counts)
    hold off
%     t = 0 : .01 : 2 * pi;

%     plot(centers, counts(:,ii));
%     xlim([-pi, pi])
%     ylim([0, max(max(counts))])
    
%     subplot(2,1,2)
%     
%     plotTime = linspace(0,112, size(circVar,1));
%     
%     plot(plotTime,circVar)
%     hold on
%     plot([plotTime(ii), plotTime(ii)], [0, 0.5], 'r')
% %     
%     hold off
%     xlim([0, 112])
%     ylim([0, 0.5])
    
    drawnow
    
end % END FOR

%% Delta phi complete

tmpDeltaPhi = squeeze(deltaPhi(:,1,:));
tmpDeltaPhi = reshape(tmpDeltaPhi', [],1);

gridPhi = squeeze(mean(deltaPhi,2));
gridPhi = reshape(gridPhi', [],1);

ax1 = subplot(2,1,1);
plot(linspace(0, 1000/60, size(gridPhi,1)), smooth(gridPhi))

gridMean = smooth(mean(p,2));

ax2 = subplot(2,1,2);
plot(linspace(0, 1000/60, size(gridMean,1)), gridMean)

linkaxes([ax1,ax2], 'x');

%% Plot VM