function [ varargout ] = GenFakeMove(  )

titleFmt = '%3.5f';

cla
xlim([0,1920])
ylim([0,1080])

hold on
plot([640, 640], [1080, 360], 'LineWidth', 10)
plot([640*2, 640*2], [1080, 360], 'LineWidth', 10)
plot([640, 640], [200, 0], 'LineWidth', 10)
plot([640*2, 640*2], [200, 0], 'LineWidth', 10)
hold off


% get current figure event functions
currFcn = get(gcf, 'windowbuttonmotionfcn');
currFcn2 = get(gcf, 'windowbuttondownfcn');
currTitle = get(get(gca, 'Title'), 'String');

% add data to figure handles
handles = guidata(gca);
if (isfield(handles,'ID') & handles.ID==1)
    disp('gtrack is already active.');
    return;
else
    handles.ID = 1;
end
handles.currFcn = currFcn;
handles.currFcn2 = currFcn2;
handles.currTitle = currTitle;
handles.theState = uisuspend(gcf);
guidata(gca, handles);

% set event functions
set(gcf,'Pointer','crosshair');
set(gcf, 'windowbuttonmotionfcn', @MouseMove);
set(gcf, 'windowbuttondownfcn', @MouseDown);    


% declare variables
xInd = 0;
yInd = 0;
mouseData = [];


    function MouseMove(src,evnt)
        
        % get mouse position
        pt = get(gca, 'CurrentPoint');
        xInd = pt(1, 1);
        yInd = pt(1, 2);
        
        % check if its within axes limits
        xLim = get(gca, 'XLim');
        yLim = get(gca, 'YLim');
        if xInd < xLim(1)
            xInd = xLim(1);
        elseif xInd > xLim(2)
            xInd = xLim(2);
        end
        
        if yInd < yLim(1)
            yInd = yLim(1);
        elseif yInd > yLim(2)
            yInd = yLim(2);
        end
        
        mouseData(end+1, :) = [xInd, yInd];
        assignin('base','mouseData',mouseData)
        
        % update figure title
        try
            title(['X = ' num2str(xInd,titleFmt) ', Y = ' num2str(yInd,titleFmt)]);
            % possibility of wrong format strings...
        catch
            error('GTRACK: Error printing coordinates. Check that you used a valid format string.')
        end
        
    end

    function MouseDown(src,evnt)
        
        % if left button, terminate
        if strcmp(get(gcf,'SelectionType'),'alt')
            varargout{1} = mouseData;
            return
        end
        
    end
end