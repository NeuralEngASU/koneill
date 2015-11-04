function [ output_args ] = GridPlot( data, layout)

chans = min(layout(layout(:)~=-1)) : max(layout(:));
border = 0.05;

mapCol = grey(128);

% Find Lower Left corner
[~,chanLLIdx] = intersect(layout,chans);

% Subdivide each section into [x,y] sections. Where x and y are layout
% dimensions.
rowsFix = linspace(0+border,1-border,size(layout,1)+1);

colsFix = linspace(0+border,1-border,size(layout,2)+1);
clf;
figure;

for ii = 1:length(chanLLIdx)
    
    chanIdx = layout(chanLLIdx(ii));
    [offx, offy] = ind2sub([size(layout,1),size(layout,2)],find(layout==chanIdx));
    
    for jj = 1:1%length(chanLLIdx)
        
        chanIdx2 = layout(chanLLIdx(jj));
        
        [xpos, ypos] = ind2sub([size(layout,1),size(layout,2)],find(layout==chanIdx2));
        
        xSubPos = [rowsFix(1,xpos), rowsFix(1,xpos+1), rowsFix(1,xpos+1), rowsFix(1,xpos)  ] + offx;
        ySubPos = [colsFix(1,ypos), colsFix(1,ypos)  , colsFix(1,ypos+1), colsFix(1,ypos+1)] + offy;
        
        colorPatch = mapCol(floor(data(chanIdx2) * 127+1),:);
        
        pData  = patch( xSubPos, ySubPos, colorPatch);
        
        
        xlim([0.5,size(layout,1)+1.5])
        ylim([0.5,size(layout,2)+1.5])
        
        set(pData, 'EdgeColor', 'none')
        
    end % END FOR channel
end % END FOR num LL corners
end % END FUNCTION

