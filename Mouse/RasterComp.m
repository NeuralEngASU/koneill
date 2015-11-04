%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	RasterComp
%		John Dylan Wright and Kevin O'Neill
%		Mario Capecchi Lab
%		20131103
%		v0.1
%		PI: Naveen Nagarajan
%
%	Inputs:
%		SpikeData1: Data loaded from the expYYYY-MM-DD_hh-mm-ss_SE.mat file
%											(nse2mat output)
%
%		SpikeData2: Data loaded from the expYYYY-MM-DD_hh-mm-ss_SE.mat file
%											(nse2mat output)
%
%       saveData: Boolean. If true output image files.
%
%       channels: The channels that will be used in the comparison
%
%       stimOn: Boolean. Determines if a backgroud patch will be used to
%               indicate when the trial/stimulation was active
%
%       stimTimes: 1x2 Vector. Needs to be altered to allow for multiple
%                  stimulation times. Should become an nx2 where the
%                  columns are time for the stimulation or trial turning on
%                  and off.
%
%       stimColor: String. Allows the user to select a color of the
%                  background patch. Uses the LineSpec color codes.
%
%		timeBounds:	Plot between the times specified (inclusive) around an event. Must be positive if events is empty.
%					Example:   [100, 1000] % seconds, no events
%					Example:   [-25, 25]   % Seconds, plot from -25 before and 25 after event
%
%		events:		If events is not empty, plot around event (such as stim turned on) and use timeBounds
%					as the plot boundaries. Events are a vector of timeStamps referencing when stimulation
%					turned on/off (give only one side of the simulation. timeBounds will tell the plotter
%					what time section around the event to plot. Units are in seconds
%					Example: trueEvents = [25, 50, 75, 100, 125, 150, ...]   on/off, on/off, on/off
%						     events = [ 25, 75, 125]    stim on events only
%							 events = [ 50, 100, 150]   stim off events only
%							 timeBounds = [-50, 100]    plot time from [event-25] to [event+100]
%
%		fs:			Sampling rate in samples/sec. Used to normalize time.
%
%       
%
%	Outputs:
%		There are no function outputs. But there is a saved file output.
%		Figures
%
%	To Use:
%		Run function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = RasterComp(SpikeData1, SpikeData2, saveData, channels, stimOn,stimTimes, stimColor,timeBounds, events, fs, varargin )

if nargin < 10
    fprintf('Not enough input arguments to MousePSTH. Quitting program. DEBUG_nargin\n');
    return
end %END IF

trialOverlay = 0.1;       % Percentage overlap between trials, try to keep it between 80% and -20%
H            = figure(1); % Initializes figure
plotColor{1} = 'b';       % Colors for the two sets of data being compared
plotColor{2} = 'r';

spikeTimes1 = SpikeData1.timeStamps; % Cell containing the spike times for all channels
spikeTimes2 = SpikeData2.timeStamps; % Cell containing the spike times for all channels

startTime1  = events(1) + timeBounds(1);	% Seconds
endTime1    = events(1) + timeBounds(2);    % Seconds

startTime2  = events(2) + timeBounds(1);    % Seconds
endTime2    = events(2) + timeBounds(2);    % Seconds

startTime1 = startTime1 * 10e5; % micro seconds
endTime1   = endTime1 * 10e5;	% micro seconds

startTime2 = startTime2 * 10e5;	% micro seconds
endTime2   = endTime2 * 10e5;	% micro seconds

% Finds spikes between the start and end times
trialSpikesIdx1 = spikeTimes1{channels(1)} > startTime1  & spikeTimes1{channels(1)} < endTime1;
trialSpikesIdx2 = spikeTimes2{channels(2)} > startTime2  & spikeTimes2{channels(2)} < endTime2;

trialSpikes1 = spikeTimes1{channels(1)}(trialSpikesIdx1);
trialSpikes2 = spikeTimes2{channels(2)}(trialSpikesIdx2);

% Normalizes the times
trialSpikes1 = trialSpikes1 - events(1)* 10e5;
trialSpikes2 = trialSpikes2 - events(2)* 10e5;

trialSpikes{1,1} = trialSpikes1;
trialSpikes{2,1} = trialSpikes2;

hold on

% Plot background patch
if stimOn
    for q = 1:length(stimTimes/2)
        bgx = [stimTimes(1)-events(1), stimTimes(2)-events(1)];
        bgxx = [bgx, fliplr(bgx)];
    end % END FOR
end % END IF

bgData = patch(bgxx, [20-10*trialOverlay, 20-10*trialOverlay, 0, 0], 1);

for i=1:length(trialSpikes)
    % Checks if there were spikes within the trial region
    if ~isempty(trialSpikes{i})
        spikeNum = length(trialSpikes{i});
        
        xx = ones(1, 3 * spikeNum) * nan;
        yy = ones(1, 3 * spikeNum) * nan;
        
        % Convert to seconds and organize xx in to sets of 3 points
        xx(2:3:3*spikeNum) = double(trialSpikes{i}).*(1/10e5);
        xx(1:3:3*spikeNum) = xx(2:3:3*spikeNum) - eps(1);
        xx(3:3:3*spikeNum) = xx(2:3:3*spikeNum) + eps(1);
        
        yy(2:3:3*spikeNum) = (i) * 10 - 10 * trialOverlay * (i-1);
        yy(1:3:3*spikeNum) = ((i-1)) * 10 - 10 * trialOverlay * (i-1);
        
        plot(xx, yy, plotColor{i})
        
    end % END IF
end % END FOR


xlim([xx(1), xx(end)])
ylim([0, (i) * 10 - 10 * trialOverlay])

set(bgData, 'FaceColor', stimColor(1))
set(bgData, 'EdgeColor', 'none')
set(bgData, 'FaceAlpha', 0.25)

set(gca, 'XTick', [0])
set(gca, 'XTickLabel',{''})

set(gca, 'YTick', [0])
set(gca, 'YTickLabel',{''})

set(gca, 'Visible', 'Off')

set(gcf, 'Units', 'inches')
set(gcf, 'Position',[0 0 3 1.5])
set(gcf, 'PaperUnits','inches','PaperPosition',[0 0 3 1.5])

shortName = ['RasterComp_', SpikeData1.expDate, 'x', SpikeData2.expDate,'_chan', num2str(channels(1)),'_chan', num2str(channels(2))];

if saveData
    print('-dpng', ['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.png'], '-r100');
    saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.fig']);
    %			saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.png']);
    saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.eps'],'epsc2');
end % END IF

end % END FUNCTION

% EOF