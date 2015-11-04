

function [] = PlotPLIEvo(pathName, g, p2, chanPairs, type, mapCol)

if isempty(type); type = 'space-invader'; end % END IF

PLI = p2';

numChans = 2^nextpow2(max(chanPairs(:)));

tmpPLI = ones(numChans,numChans) * -1;
tmpPLI(sub2ind([numChans,numChans], chanPairs(:,1), chanPairs(:,2))) = PLI;
tmpPLI(sub2ind([numChans,numChans], chanPairs(:,2), chanPairs(:,1))) = PLI;
tmpPLI(logical(eye(numChans))) = zeros(1,numChans);

colorMat = mapCol;
colormap(mapCol);

if strcmp(type, 'confusion')
    
    imagesc(tmpPLI, [0, 0.7]);
    
elseif strcmp(type, 'space-invader')
    %%
    layout = g.layout;
    
    chans = 1:numChans;
    border = 0.05;
    
    % Find Lower Left corner
    [~,chanLLIdx] = intersect(layout,chans);
        
    % Subdivide each section into [x,y] sections. Where x and y are layout
    % dimensions.
    rowsFix = linspace(0+border,1-border,size(layout,1)+1);
    
    colsFix = linspace(0+border,1-border,size(layout,2)+1);
    clf;
    figure(1);
    
    for ii = 1:length(chanLLIdx)
        
        chanIdx = layout(chanLLIdx(ii));
        [offx, offy] = ind2sub([size(layout,1),size(layout,2)],find(layout==chanIdx));
                
        if isempty(intersect(chanIdx, g.badchan))
            
            for jj = 1:length(chanLLIdx)
                
                chanIdx2 = layout(chanLLIdx(jj));
                
                [xpos, ypos] = ind2sub([size(layout,1),size(layout,2)],find(layout==chanIdx2));
                
                xSubPos = [rowsFix(1,xpos), rowsFix(1,xpos+1), rowsFix(1,xpos+1), rowsFix(1,xpos)  ] + offx;
                ySubPos = [colsFix(1,ypos), colsFix(1,ypos)  , colsFix(1,ypos+1), colsFix(1,ypos+1)] + offy;
                
%                 if tmpPLI(chanIdx, chanIdx2) > 0
%                     colorPatch = colorMat(end,:);
%                 elseif tmpPLI(chanIdx, chanIdx2) < 0
%                     colorPatch = colorMat(1,:);
%                 else
%                     colorPatch = colorMat(128,:);
%                 end
%                 
                colorPatch = colorMat(floor(tmpPLI(chanIdx,chanIdx2) * 127+1),:);
                
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
    
    title(sprintf('%s: PLI', g.subject))
    
    xlim([0.9,size(layout,1)+1.1])
    ylim([0.9,size(layout,2)+1.1])
    
    axis('square')
    camroll(-90)
    set(gca, 'XTick', [])
    set(gca, 'XTickLabel', [])
    set(gca, 'YTick', [])
    set(gca, 'YTickLabel', [])
        
    cData = colorbar('Location', 'East');
    
    set(cData, 'YAxisLocation','right')
    set(cData, 'YTick', [0, 0.25,0.5,0.75, 1])
    set(cData, 'YTickLabel', [0, 0.25,0.5,0.75, 1])
    set(cData, 'TickDirection', 'Out')
%     set(cData, 'Label', 'PLI')
    cData.Label.String = 'PLI';
    
%     savefig(fullfile(pathName, sprintf('%s_Sig', g.subject)))
%     print(fullfile(pathName, sprintf('%s_Sig', g.subject)), '-dpng')
    
    
end % END IF

end % END FUNCTION

% EOF