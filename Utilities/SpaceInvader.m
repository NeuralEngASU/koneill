function [  ] = SpaceInvader( outputPath, outputName, gridDef, data, varargin )

%% Parse Inputs
% Defaults
LINEARFLAG = 0; % If the data formated in a linear form
CHANPAIRS = []; % Pairs of channels
PLOTDATA = struct('border', 0.05,...
    'title', '',...
    'xlabel', '',...
    'ylabel', '',...
    'xtick', [],...
    'ytick', [],...
    'xtickdir', 'in',...
    'ytickdir', 'in',...
    'ticklength', [0,0],...
    'cbarlabel', '',...
    'cbartick', [0:0.25:1],...
    'cbarticklabel', [0:0.25:1],...
    'cbartickdir', 'none');
DESIREDGRID = 1; % Desired grid to plot
MINVAL = 0;
MAXVAL = 1;

MAPCOL = jet(128);

% Parse Varargin
for ii = 1:2:length(varargin)
    if ~exist(upper(varargin{ii}), 'var')
        fprintf('Unknown option entry: %s\n', varargin{ii})
        return;
    else
        eval([upper(varargin{ii}) '=varargin{ii+1};']);
    end % END IF exist
end % END FOR varargin

linearFlag = LINEARFLAG;
chanPairs = CHANPAIRS;
plotData = PLOTDATA;
desiredGrid = DESIREDGRID;
mapCol = MAPCOL;
minVal = MINVAL;
maxVal = MAXVAL;

clear LINEARFLAG CHANPAIRS PLOTDATA DESIREDGRID MAPCOL

% If chanPairs was not given
if linearFlag && isempty(chanPairs) && ~isfield(gridDef, 'chanPairs')
    chanPairs = nchoosek(1:gridDef.numChan,2);
end % END IF
%%

minVal = min(data(:));
maxVal = max(data(:));

data = (data+minVal);
data = data./(max(data(:)));

% data = data(:);

if linearFlag
    
    data = data(:);
    
    numChans = gridDef.numChan;
    
    tmpData = ones(numChans,numChans) * -1;
    tmpData(sub2ind([numChans,numChans], chanPairs(:,1), chanPairs(:,2))) = data;
    tmpData(sub2ind([numChans,numChans], chanPairs(:,2), chanPairs(:,1))) = data;
    tmpData(logical(eye(numChans))) = zeros(1,numChans);
    
    data = tmpData;
    clear tmpData
    
end % END IF linearFlag

%% Make Plot

desiredGrid = desiredGrid(:);

for kk = 1:length(desiredGrid(:))
    gridNum = desiredGrid(kk);
    
    if iscell(gridDef.layout)
        layout = gridDef.layout{gridNum};
    else
        layout = gridDef.layout;
    end % END IF iscell
    
    chans = min(layout(layout(:)~=-1)) : max(layout(:));
    border = 0.05;
    
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
        
        if isempty(intersect(chanIdx, gridDef.badChan))
            
            for jj = 1:length(chanLLIdx)
                
                chanIdx2 = layout(chanLLIdx(jj));
                
                [xpos, ypos] = ind2sub([size(layout,1),size(layout,2)],find(layout==chanIdx2));
                
                xSubPos = [rowsFix(1,xpos), rowsFix(1,xpos+1), rowsFix(1,xpos+1), rowsFix(1,xpos)  ] + offx;
                ySubPos = [colsFix(1,ypos), colsFix(1,ypos)  , colsFix(1,ypos+1), colsFix(1,ypos+1)] + offy;
                
                colorPatch = mapCol(floor(data(chanIdx,chanIdx2) * 127+1),:);
                
                pData  = patch( xSubPos, ySubPos, colorPatch);
                
                if  chanIdx == chanIdx2 || ~isempty(intersect(chanIdx2, gridDef.badChan))
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
    
    title(sprintf('%s', gridDef.subject))
    
    xlim([0.9,size(layout,1)+1.1])
    ylim([0.9,size(layout,2)+1.1])
    
    axis('square')
    camroll(-90)
    set(gca, 'XTick', [])
    set(gca, 'XTickLabel', [])
    set(gca, 'YTick', [])
    set(gca, 'YTickLabel', [])
    colormap(jet)
    cData = colorbar('Location', 'East');
    
    set(cData, 'YAxisLocation','right')
    set(cData, 'YTick', plotData.cbartick)
    set(cData, 'YTickLabel', plotData.cbarticklabel)
    set(cData, 'TickDirection', 'Out')
    %     set(cData, 'Label', 'PLI')
    cData.Label.String = plotData.cbarlabel;
    
end % END FOR gridNum
end % END FUNCTION

% EOF