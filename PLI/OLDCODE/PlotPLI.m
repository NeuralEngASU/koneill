

function [] = PlotPLI(path, g, p, chanPairs, type)

if isempty(type); type = 'confusion'; end % END IF

PLI = squeeze(mean(p,1));

PLI = PLI(:,1);

numChans = 2^nextpow2(max(chanPairs(:)));

tmpPLI = ones(numChans,numChans) * -1;
tmpPLI(sub2ind([numChans,numChans], chanPairs(:,1), chanPairs(:,2))) = PLI;
tmpPLI(sub2ind([numChans,numChans], chanPairs(:,2), chanPairs(:,1))) = PLI;
tmpPLI(logical(eye(numChans))) = zeros(1,numChans);

colorMat = jet(128);

if strcomp(type, 'confusion')
    
    imagesc(tmpPLI, [0, 0.15]);
    
elseif strcomp(type, 'space-invader')
    %%
    layout = g.layout;
    
    chans = 1:numChans;
    border = 0.05;
    
    % Find Lower Left corner
    [~,chanLLIdx] = intersect(layout,chans);
    
    [rows,cols] = ind2sub([size(layout,1),size(layout,2)], chanLLIdx);
    
    % Subdivide each section into [x,y] sections. Where x and y are layout
    % dimensions.
    
%     [rowsFix, meshRows] = meshgrid(linspace(0+border,1-border,size(layout,1)+1), rows);
%     rows = meshRows + rowsFix;
%     
%     [colsFix, meshCols] = meshgrid(linspace(0+border,1-border,size(layout,2)+1), cols);
%     cols = meshCols + colsFix;
    
    rowsFix = linspace(0+border,1-border,size(layout,1)+1);
    
    colsFix = linspace(0+border,1-border,size(layout,2)+1);
    
    for ii = 1:length(chanLLIdx)
        
        chanIdx = layout(chanLLIdx(ii));
        [offx, offy] = ind2sub([size(layout,1),size(layout,2)],find(layout==chanIdx));
                
        if isempty(intersect(chanIdx, g.badchan))
            
            for jj = 1:length(chanLLIdx)
                
                chanIdx2 = layout(chanLLIdx(jj));
                
                [xpos, ypos] = ind2sub([size(layout,1),size(layout,2)],find(layout==chanIdx2));
                
                xSubPos = [rowsFix(1,xpos), rowsFix(1,xpos+1), rowsFix(1,xpos+1), rowsFix(1,xpos)  ] + offx;
                ySubPos = [colsFix(1,ypos), colsFix(1,ypos)  , colsFix(1,ypos+1), colsFix(1,ypos+1)] + offy;
                
                colorPatch = colorMat(floor(tmpPLI(chanIdx,chanIdx2)/0.15 * 127+1),:);
                
                pData  = patch( xSubPos, ySubPos, colorPatch);
                
                if  chanIdx == chanIdx2 || ~isempty(intersect(chanIdx2, g.badchan))
                    set(pData, 'FaceColor', [0,0,0])
                end % END IF
                
                xlim([0.5,size(layout,1)+1.5])
                ylim([0.5,size(layout,2)+1.5])
                
                set(pData, 'EdgeColor', 'none')
                
            end % END FOR
            
            blackOut = find(layout == -1);
            
            for blackIdx = 1:length(blackOut)
            
                [xpos, ypos] = ind2sub([size(layout,1),size(layout,2)],blackOut(blackIdx));
                
                xSubPos = [rowsFix(1,xpos), rowsFix(1,xpos+1), rowsFix(1,xpos+1), rowsFix(1,xpos)  ] + offx;
                ySubPos = [colsFix(1,ypos), colsFix(1,ypos)  , colsFix(1,ypos+1), colsFix(1,ypos+1)] + offy;
                                
                pData  = patch( xSubPos, ySubPos, colorPatch);
                set(pData, 'FaceColor', [0,0,0])
                
                xlim([0.5,size(layout,1)+1.5])
                ylim([0.5,size(layout,2)+1.5])
                
                set(pData, 'EdgeColor', 'none')
                
            end % END FOR blackIdx
            
            
            [offx, offy] = ind2sub([size(layout,1),size(layout,2)],find(layout==chanIdx));
            
            pData = patch([rowsFix(1), rowsFix(end), rowsFix(end),   rowsFix(1)] + offx,...
                          [colsFix(1),   colsFix(1), colsFix(end), colsFix(end)] + offy,...
                [1,1,1]);
            
            set(pData, 'FaceColor', 'none')
            set(pData, 'EdgeColor', [0,0,0]);
            
        end % END IF
        
    end % END FOR
    
    axis('square')
    camroll(-90)
    set(gca, 'XTick', [])
    set(gca, 'XTickLabel', [])
    set(gca, 'YTick', [])
    set(gca, 'YTickLabel', [])
    
    
end % END IF

end % END FUNCTION

% EOF