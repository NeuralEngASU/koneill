% Mono-Phasic Data Plot


%% TODO
    % Plot
    % Constrain y axis zoom
    % Background Patches
    % Vertical Deliminators
    % Scaling

%% Test Plot

numChans = 16;

data = ones(numChans, 1000).*repmat(sind(2*pi*(1:1000)), numChans,1);

offset = repmat([10:10:numChans*10]', 1, 1000);

dataOff = data + offset;

%%
h = plot(1:1000,dataOff');

xlim([0, 1000])
ylim([0, numChans*10+10])

set(gca, 'YTick', [10:10:numChans*10])
set(gca, 'YTickLabel', [1:numChans])

% EOF