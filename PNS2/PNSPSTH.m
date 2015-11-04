%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	MousePSTH
%		Kevin O'Neill
%		Mario Capecchi Lab
%		20131103
%		v0.1
%		PI: Naveen Nagarajan
%
%	Inputs:
%		SpikeData: 	Data loaded from the expYYYY-MM-DD_hh-mm-ss_SE.mat file
%											(nse2mat output)
%
%		saveDate:  	Saves plots or regenerated plots to the hard drive. Data will be saved as
%					a .fig, .png, and a .eps. File names will follow the convention:
%						PSTH_YYYY-MM-DD_hh-mm-ss.ext
%					saveData = 1: Save data
%					saveData = 0: Do not save data
%
%		stimOn:		Plots an overlay of the stimulation as a background color (blue and yellow)
%					stimOn = 1: plot overlays
%					stimOn = 0: don't plot overlays (Default)
%
%		stimTimes:	Holds the time stamps for when the stimulus turned on and off. Must contain 2 entries
%					for stim on and stim off.
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
%	Outputs:
%		There are no function outputs. But there is a saved file output.
%		Figures
%
%	To Use:
%		Run function. All inputs must be given.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = MousePSTH( SpikeData, saveData, stimOn, stimTimes, stimColor, timeBounds, events, fs)

if nargin < 8
    fprintf('Not enough input arguments to MousePSTH. Quitting program. DEBUG_nargin\n');
    return
end %END IF

buffer = 2*fs; % Buffer for the beginning and end data in samples. Used to avoid edge effects.

spikeTimes = SpikeData.timeStamps; % Cell containing the spike times for all channels

if timeBounds(1) < 2
    fprintf('timeBounds must be at least 2 seconds away from the beginning and end of the file. DEBUG_time\n');
end % END IF

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

numChans = size(spikeTimes,1);

for j = 1:numChans
    for i = 1:numEvents
        
        startTime = events(i) + timeBounds(1);	% Seconds
        endTime   = events(i) + timeBounds(2);	% Seconds
        
        % 			startTime = startTime * fs;				% Samples
        % 			endTime   = endTime * fs;				% Samples
        
        startTime = startTime * 10e5;			% micro seconds
        endTime   = endTime * 10e5;				% micro seconds
        
        % Finds the indices of spikes that occur between [startTime - buffer] and [endTime + buffer]
        trialSpikesIdx = spikeTimes{j} > startTime - buffer & spikeTimes{j} < endTime + buffer;
        
        % Selects the timeStamps using the indices found above
        trialSpikes = spikeTimes{j}(trialSpikesIdx);
        
        % Normalizes the time to 0 seconds
        trialSpikes = trialSpikes - startTime;
        
        % Spike times need to be placed into a structure for psth() to function.
        AP.spikeTimes = double(trialSpikes)./10e5;
        
        %Counts how many spikes are in the selected time slice
        numSpikes = length(trialSpikes);
        
        % Ensures that only the channels with spikes are analysed
        if numSpikes
            
            figure(1)
            % Preforms the PSTH (from Chronux)
            [R, t, E] = psth(AP, 0.300, 'r', [0,3], 1);
            
            close(gcf)
            % Makes sure all vectors are vertical
            sz = size(R);
            if sz(1) < sz(2)
                R=R'; t=t'; E=E';
            end % END IF
            
            % Stores trials into cells
            RCell{i} = R;	% mean
            tCell{i} = t;	% time
            ECell{i} = E;	% error
            
        else
            RCell{i} = NaN;
            tCell{i} = NaN;
            ECell{i} = NaN;
        end % END IF
        
        
    end % END FOR
    
    % Removes NaN information
    count = 1;
    for c = 1:length(RCell)
        
        if ~isnan(RCell{c})
            RMat(:,count) = RCell{c};
            EMat(:,count) = ECell{c};
            count = count+1;
        end
        
    end % END FOR
    
    % Averages the firing rate for all events
    RMat = mean(RMat, 2);
    EMat = mean(EMat, 2);
    
    % Stores the firing rate and error for all channels
    if j == 1
        rateMat = zeros(length(RMat), numChans, 2);
    end;
    rateMat(:, j, 1) = RMat;
    rateMat(:, j, 2) = EMat;
    
end % END FOR

%% PSTH Plot

% Creates a time vector around the event (typically around 0)
x = linspace(timeBounds(1)-2, timeBounds(2)+2, length(rateMat(:,1,1)));
xx = [x, fliplr(x)];

if stimOn
    for q = 1:length(stimTimes/2)
        bgx = [stimTimes(1), stimTimes(2)];
        bgxx = [bgx, fliplr(bgx)];
    end % END FOR
end % END IF

for k = 1:size(rateMat,2)
    
    
    dataTop = rateMat(:,k,1)+rateMat(:,k,2);
    dataBot = rateMat(:,k,1)-rateMat(:,k,2);
    dataMid = rateMat(:,k,1);
    
    data1 = [[rateMat(:,k,1)+rateMat(:,k,2)]', fliplr([rateMat(:,k,1)-rateMat(:,k,2)]')];
    
    figure(k)
    hold on
    lData  = plot(x, dataMid, '-b');
    bgData = patch(bgxx, [2*max(data1), 2*max(data1), 0, 0], 1);
    pData  = patch(xx, data1, 1);
    lData  = plot(x, dataMid, '-b');
    hold off
    
    set(bgData, 'FaceColor', stimColor(1))
    set(bgData, 'EdgeColor', 'none')
    set(bgData, 'FaceAlpha', 0.25)
    
    set(pData, 'FaceColor', 'k')
    set(pData, 'EdgeColor', 'none')
    set(pData, 'FaceAlpha', 0.25)
    
    set(lData, 'LineWidth', 2)
    
    xlim([timeBounds(1), timeBounds(2)])
    ylim([0, 1.1*max(data1)])						%***** Hard Code ylim values if needed *****%
    
    title(sprintf('Instantaneous Firing Rate\nChannel: %d', k))
    xlabel('Time, sec')
    ylabel('AP Firing Rate, Hz')
    
    set(gcf, 'Units', 'inches')
    set(gcf, 'Position',[0 0 5 3])
    set(gcf, 'PaperUnits','inches','PaperPosition',[2 2 5 3])
    
    legend('Firing Rate', 'Location', [0.89, 0.8, 0.125, 0.0625])
    
    shortName = ['PSTH_', SpikeData.expDate, '_chan', num2str(k)];
    
    if saveData
        print('-dpng', ['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.png'], '-r100');
        saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.fig']);
        %			saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.png']);
%         saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.eps'],'epsc2');
        saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.pdf']);

    end % END IF
    
end % END FOR

end % END FUNCTION


% EOF