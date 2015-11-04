

function [] = PlotPLI2(pathName, g, p, r, chanPairs, type, maxPLI, maxR)

if isempty(type); type = 'space-invader'; end % END IF

PLI = squeeze(mean(p,1));

PLI = PLI(:,1);

numChans = 2^nextpow2(max(chanPairs(:)));

tmpPLI = ones(numChans,numChans) * -1;
tmpPLI(sub2ind([numChans,numChans], chanPairs(:,1), chanPairs(:,2))) = PLI;
tmpPLI(sub2ind([numChans,numChans], chanPairs(:,2), chanPairs(:,1))) = PLI;
tmpPLI(logical(eye(numChans))) = zeros(1,numChans);

colorMat = jet(128);

R = squeeze(mean(r,1));

R = R(:,1);

tmpR = ones(numChans,numChans) * -1;
tmpR(sub2ind([numChans,numChans], chanPairs(:,1), chanPairs(:,2))) = R;
tmpR(sub2ind([numChans,numChans], chanPairs(:,2), chanPairs(:,1))) = R;
tmpR(logical(eye(numChans))) = ones(1,numChans);

<<<<<<< HEAD

colormap('jet')

if strcmp(type, 'confusion')
    
    imagesc(tmpPLI, [0, 0.15]);
    
elseif strcmp(type, 'space-invader')
=======
if strcomp(type, 'confusion')
    
    imagesc(tmpPLI, [0, 0.15]);
    
elseif strcomp(type, 'space-invader')
>>>>>>> origin/master
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
                
                colorPatch = colorMat(floor(tmpPLI(chanIdx,chanIdx2)/maxPLI * 127+1),:);
                
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
<<<<<<< HEAD
        
=======
    
>>>>>>> origin/master
    cData = colorbar('Location', 'East');
    
    set(cData, 'YAxisLocation','right')
    set(cData, 'YTick', [0, 1])
    set(cData, 'YTickLabel', [0, maxPLI])
    set(cData, 'TickDirection', 'Out')
<<<<<<< HEAD
%     set(cData, 'Label', 'PLI')
    cData.Label.String = 'PLI';
    
    savefig(fullfile(pathName, sprintf('%s_PLI', g.subject)))
    print(fullfile(pathName, sprintf('%s_PLI', g.subject)), '-dpng')
=======
    set(cData, 'Label', 'PLI')
%     set(cData, 'YLabel', 'PLI')

    
    save(fullfile(pathName, sprintf('%s_PLI', g.subject)))
    print(fullfile(pathName, sprintf('%s_PLI', g.subject)), '-png')
>>>>>>> origin/master
    
    
    %%
    figure(2)
    
    for ii = 1:length(chanLLIdx)
            
        chanIdx = layout(chanLLIdx(ii));
        [offx, offy] = ind2sub([size(layout,1),size(layout,2)],find(layout==chanIdx));
                
        if isempty(intersect(chanIdx, g.badchan))
            
            for jj = 1:length(chanLLIdx)
                
                chanIdx2 = layout(chanLLIdx(jj));
                
                [xpos, ypos] = ind2sub([size(layout,1),size(layout,2)],find(layout==chanIdx2));
                
                xSubPos = [rowsFix(1,xpos), rowsFix(1,xpos+1), rowsFix(1,xpos+1), rowsFix(1,xpos)  ] + offx;
                ySubPos = [colsFix(1,ypos), colsFix(1,ypos)  , colsFix(1,ypos+1), colsFix(1,ypos+1)] + offy;
                
                colorPatch = colorMat(floor(tmpR(chanIdx,chanIdx2)/maxR * 127+1),:);
                
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
    
    title(sprintf('%s: R', g.subject))
    
    xlim([0.9,size(layout,1)+1.1])
    ylim([0.9,size(layout,2)+1.1])
    
    axis('square')
    camroll(-90)
    set(gca, 'XTick', [])
    set(gca, 'XTickLabel', [])
    set(gca, 'YTick', [])
    set(gca, 'YTickLabel', [])
<<<<<<< HEAD
        
    cData = colorbar('Location', 'East');
    
    set(cData, 'YAxisLocation','right')
    set(cData, 'YTick', [0, 1])
    set(cData, 'YTickLabel', [0, maxR])
    set(cData, 'TickDirection', 'Out')
%     set(cData, 'Label', 'PLI')
    cData.Label.String = 'R';
    
    savefig(fullfile(pathName, sprintf('%s_R', g.subject)))
    print(fullfile(pathName, sprintf('%s_R', g.subject)), '-dpng')
=======
    
    save(fullfile(pathName, sprintf('%s_R', g.subject)))
    print(fullfile(pathName, sprintf('%s_R', g.subject)), '-png')
>>>>>>> origin/master
    
    %%    
    
end % END IF

end % END FUNCTION

% EOF