%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	TrialRaw
%		Kevin O'Neill
%		Mario Capecchi Lab
%		20131103
%		v0.1
%		PI: Naveen Nagarajan
%
%	Inputs:
%		CSCData: Data loaded from the expYYYY-MM-DD_hh-mm-ss_CSC.mat file
%											(nsc2mat output)
%
%       saveDate:  	Saves plots or regenerated plots to the hard drive. Data will be saved as
%					a .fig, .png, and a .pdf. File names will follow the convention:
%						PSTH_YYYY-MM-DD_hh-mm-ss.ext
%					saveData = 1: Save data
%					saveData = 0: Do not save data
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
%       ADDIDTIONAL NOTE: This function is currently set up to run with
%       only one channel and event input. It can be modified to loop over
%       multiple channels if needed; just follow struture from other
%       functions.
%
%	Outputs:
%		There are no function outputs. But there is a saved file output.
%		Figures
%
%	To Use:
%		Run function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = TrialRaw( CSCData, saveData, channel, stimOn, stimTimes, stimColor, timeBounds, events)

if nargin < 8
    fprintf('Not enough input arguments to TrialRaw. Quitting program. DEBUG_nargin\n');
    return
end %END IF


% Ensures that there is at least one arbitrary event
if ~isempty(events)
    numEvents = size(events);
else
    numEvents = 1;
    events = 0;
end % END IF

% If stimOn is selected, makes sure that the information is correct.
if stimOn
    if length(stimTimes) < 2
        fprintf('stimTimes not initialized correctly. Wrong size. DEBUG_stimTimes\n');
    end % END IF
    
    if length(stimColor) ~= length(stimTimes)/2
        fprintf('stimColor not initialized correctly. Wrong size. DEBUG_stimColor\n');
    end % END IF
end % END IF

timeStamps = CSCData.timeStamps;

startTime = events(1)   + timeBounds(1);	% Seconds
endTime   = events(end) + timeBounds(2);	% Seconds

startTime = startTime * 10e5;			% micro seconds
endTime   = endTime   * 10e5;			% micro seconds

timeIdx = find(timeStamps >= startTime & timeStamps <= endTime);

% Selects the timeStamps using the indices found above
tempChanName = ['CSC',num2str(channel)];
voltData = h5read(CSCData.hdf5FileName, ['/', tempChanName], [1, timeIdx(1)], [1, length(timeIdx)*512]);
voltData = voltData.*1000; % Convert to mV

time = linspace(startTime*(1/10e5), endTime*(1/10e5), length(voltData));

hold on

if stimOn
    for q = 1:length(stimTimes/2)
        bgx = [stimTimes(1)-events(1), stimTimes(2)-events(1)];
        bgxx = [bgx, fliplr(bgx)];
    end % END FOR
    bgData = patch(bgxx, [2*max(voltData), 2*max(voltData), 2*min(voltData), 2*min(voltData)], 1); 
    % Assumes the voltage crosses 0. Needs another method if the voltage is
    % one one side of the zero line.
    
    set(bgData, 'FaceColor', stimColor(1))
    set(bgData, 'EdgeColor', 'none')
    set(bgData, 'FaceAlpha', 0.25)
end % END IF

plot(time, voltData, 'k', 'linewidth', 1)

hold off

title('Raw Voltage Trace')
xlabel('Time, sec')
ylabel('Voltage, mV')

xlim([timeBounds(1), timeBounds(2)])
ylim([1.2*min(voltData), 1.2*max(voltData)])


% 
% set(gca, 'XTick', [0])
% set(gca, 'XTickLabel',{''})
% 
% set(gca, 'YTick', [0])
% set(gca, 'YTickLabel',{''})

set(gca, 'Visible', 'On')

set(gcf, 'Units', 'inches')
set(gcf, 'Position',[0 0 3 1.5])
set(gcf, 'PaperUnits','inches','PaperPosition',[0 0 3 1.5])


if stimOn
    shortName = ['RawTraceStim_', CSCData.expDate, '_chan', num2str(channel)];
else
    shortName = ['RawTrace_', CSCData.expDate, '_chan', num2str(channel)];
end % END IF

% if saveData
%     print('-dpng', ['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.png'], '-r100');
%     saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.fig']);
%     %			saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.png']);
%     saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.eps'],'epsc2');
% end % END IF

end % END FUNCTION