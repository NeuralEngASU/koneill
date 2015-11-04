%% Convolution Visualization

boxWidth = 1000;
box = [zeros(1,0), ones(1,boxWidth),zeros(1,0)];

triWidth = 1000;
tri = [zeros(1,500), -1/triWidth*(1:triWidth) + 1,zeros(1,500)];


tBox = 1:length(box);

tTri = 1:length(tri);

for ii = 1:length(tBox)*4

clf
    
plot(tBox - length(tBox) * 2 + ii, box, 'b', 'LineWidth', 2.75);
hold on
plot(tTri, tri, 'r', 'LineWidth', 2.75);

xlim([-250, 250])

drawnow
% pause(100)
end % END FOR


%%

plot(conv(box,tri), 'k', 'LineWidth', 2.75); xlim([0, 6000]); ylim([-0.1, 510]);
set(gca, 'XTick', [])
set(gca, 'YTick', [])

set(gcf, 'Units', 'inches')
set(gcf, 'Position',[2 2 5 2])
set(gcf, 'PaperUnits','inches','PaperPosition',[2 2 5 3])
