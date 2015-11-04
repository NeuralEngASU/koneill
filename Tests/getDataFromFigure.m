function getDataFromFigure()
% Plots graph and sets up a custom data tip update function
fig = figure;
a = -16; t = 0:60;
plot(t,sin(a*t))

% variable to store data points
myData = [];

% enable data cursor mode
dcm_obj = datacursormode(fig);
set(dcm_obj,'UpdateFcn',@myUpdateFcn)
set(dcm_obj, 'enable', 'on')

% do disable data cursor mode use
% set(dcm_obj, 'enable', 'off')


    function txt = myUpdateFcn(dummy, event_obj)
        % Customizes text of data tips

        % read out data point
        pos = get(event_obj,'Position');

        % store data point
        myData(end+1,:) = pos;

        % no data shown on figure
        txt = {''};

        % or
        % data also shown on figure:
        % txt = {['Time: ',num2str(pos(1))],...
        %         ['Amplitude: ',num2str(pos(2))]};
    end
end