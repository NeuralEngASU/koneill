%% Lines

x = 0:0.01:10;
y = real((-2).^x);

% %%
stdErr = 100*exp(x./10);

figure;
plot(stdErr)
% %%
lData = plot(x,y); % Plots y versus x

set(gcf, 'Units', 'Inches')
set(gcf, 'Position',[2 2 3 1.5])

% set(lData, 'Color', 'g');       % Changes the color of the line
% set(lData, 'LineStyle', '--');  % Changes the style of the line
% set(lData, 'LineWidth', 2.75);  % Adjusts the line thickness

%% Patches
figure;
stdErr = 75*exp(x./10);

xx = [x, fliplr(x)];    % x-values for each of the y values
yy = [[y + 2*stdErr], fliplr([y - 2*stdErr])]; % y-values must form a polygon


hold on
backData = patch([2, 6, 6, 2], [2000, 2000, -1000, -1000], 1); % Background patch
pData = patch(xx, yy, 1); % Plot patch
lData = plot(x,y);        % Plot line over patch
hold off


% hold on
% pData = patch(xx, yy, 1); % Plot patch
% lData = plot(x,y);        % Plot line over patch
% hold off

set(backData, 'FaceColor', 'k');   % Changes the color of the patch
set(backData, 'EdgeColor', 'none');% Changes the patch edge color
set(backData, 'FaceAlpha', 0.25);  % Changes the patch transparency

set(pData, 'FaceColor', 'b');   % Changes the color of the patch
set(pData, 'EdgeColor', 'b');   % Changes the patch edge color
set(pData, 'FaceAlpha', 0.25);  % Changes the patch transparency

set(lData, 'Color', 'b');       % Changes the color of the line
set(lData, 'LineStyle', '-');   % Changes the style of the line
set(lData, 'LineWidth', 2.75);  % Adjusts the line thickness

set(gcf, 'Units', 'Inches')
set(gcf, 'Position',[2 2 3 1.5])


%% Plot Format

title(sprintf('Electrode: %d', electrodes(k)), 'FontSize', 15)

ylabel('Firing Rate, Hz', 'FontSize', 15)
xlabel('Time, sec', 'FontSize', 15)

set(gcf, 'Units', 'Inches')
set(gcf, 'Position',[2 2 6 3])
% set(gcf, 'PaperUnits', 'FontSize', 'PaperPosition',[2 2 6 3])

aPos = get(gca, 'Position');

xlim([0, 3])
ylim([0, 40])

p1 = get(gca, 'position');
legend({'Hello', 'World'}, 'Location', [0.79, 0.8, 0.125, 0.0625])
legend({'Hello', 'World'}, 'Location', 'NorthEast')
legend boxoff

set(gca, 'yTick', [0, 10, 20, 30, 40])
set(gca, 'yTickLabel', [0, 10, 20, 30, 40])
set(gca, 'xTick', [0, 1, 2, 3])
set(gca, 'xTickLabel', [0, 1, 2, 3])

grid on
set(gca,'XGrid','off','YGrid','on','ZGrid','off')

set(gca,'FontSize',15)

h_xlabel = get(gca,'xLabel');
set(h_xlabel,'FontSize',15);

h_ylabel = get(gca,'yLabel');
set(h_xlabel,'FontSize',15);

% print('-dpng', ['C:\CodeRepo\Lab\PNS\P201301\Figures2\PSTH_V2_', num2str(k), '.png'], '-r100');
% saveas(gcf,['C:\CodeRepo\Lab\PNS\P201301\Figures2\PSTH_V2_', num2str(k), '.png']);




