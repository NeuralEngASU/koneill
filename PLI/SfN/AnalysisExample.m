
%%
chanIdx = 35;

% Find Channels to Plot (3x3 grid around chanIdx)
gridDim = [8,8];
[ii, jj] = ind2sub(gridDim, chanIdx);

iiNew = [ii-1, ii, ii+1];
jjNew = [jj-1, jj, jj+1];

% Find index pairs
[idx2, idx1] = find(true(numel(iiNew),numel(jjNew))); 
indPairs = [reshape(iiNew(idx1), [], 1), reshape(jjNew(idx2), [], 1)];
    
% Remove invalid pairs (where there is a 0 index)
indParisIdx = (indPairs(:,1) > 0) & (indPairs(:,2) > 0) & (indPairs(:,1) <= gridDim(1)) & (indPairs(:,2) <= gridDim(2));
indPairs2 = indPairs(indParisIdx,:);

% Convert index notation to channel number
chansPlot = sort(sub2ind([8,8], indPairs2(:,1), indPairs2(:,2)), 'ascend');

[chansPlotOrig(1), chansPlotOrig(2)] = ind2sub(gridDim, chanIdx);
chansPlot = chansPlot(chansPlot~=chanIdx);

chansPlotCoord = [];
[chansPlotCoord(:,1), chansPlotCoord(:,2)] = ind2sub(gridDim, chansPlot);



%%


figure;
hold on
mapCol = gray(128);
border = 0;

layout = reshape(1:64,8,8);

for ii = 1:4
       
    % Find Lower Left corner
    chans = 1:64;
    [~,chanLLIdx] = intersect(layout,chans);
    
    % Subdivide each section into [x,y] sections. Where x and y are layout
    % dimensions.
    rowsFix = linspace(0+border,1-border,size(layout,1)+1);
    
    colsFix = linspace(0+border,1-border,size(layout,2)+1);
    
    for jj = 1:length(chanLLIdx)
        
        chanIdx = layout(chanLLIdx(jj));
        [offx, offy] = ind2sub([size(layout,1),size(layout,2)],find(layout==chanIdx));
        
        xSubPos = [rowsFix(1,1), rowsFix(1,end), rowsFix(1,end), rowsFix(1,1)  ] + offx;
        ySubPos = [colsFix(1,1), colsFix(1,1)  , colsFix(1,end), colsFix(1,end)] + offy;
        
        chanIdx2 = layout(chanLLIdx(jj));
        
        [xpos, ypos] = ind2sub([size(layout,1),size(layout,2)],find(layout==chanIdx2));
        
        colorPatch = mapCol(floor(128),:);
        
        pData  = patch( xSubPos, ySubPos, colorPatch);
        
        
        xlim([1,size(layout,1)+1])
        ylim([1,size(layout,2)+1])
        
        set(pData, 'EdgeColor', 'k')
        axis square
                
        set(gca, 'xtick', [1.5:8.5])
        set(gca, 'xticklabel', [1:8])
        
        set(gca, 'ytick', [1.5:8.5])
        set(gca, 'yticklabel', [1:8:57])
        
        xlabel('Channel')
        ylabel('Channel')
        box on
        
    end % END FOR num LL corners    
end

% figure;

for kk = 1:size(chansPlotCoord,1)
    
    xValLine = [chansPlotCoord(kk,1), chansPlotOrig(1)] + 0.5;
    yValLine = [chansPlotCoord(kk,2), chansPlotOrig(2)] + 0.5;
    
    xValMark = chansPlotCoord(kk,1) + 0.5;
    yValMark = chansPlotCoord(kk,2) + 0.5;
    
    plot(xValLine,yValLine,'-r', 'linewidth', 2)
    plot(xValMark,yValMark,'ko', 'linewidth', 1, 'markersize', 10, 'MarkerEdgeColor','k', 'markerfacecolor', [0.8500    0.3250    0.0980])
       
end

plot(chansPlotOrig(1)+0.5, chansPlotOrig(2)+0.5, 'bo', 'linewidth', 1, 'markersize', 15, 'MarkerEdgeColor','k', 'markerfacecolor', [0    0.4470    0.7410])

