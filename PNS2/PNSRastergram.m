%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Rastergram
%		Kevin O'Neill
%		Mario Capecchi Lab
%		20131103
%		v0.1
%		PI: Naveen Nagarajan
%
%	Inputs:
%		SpikeData: Data loaded from the expYYYY-MM-DD_hh-mm-ss_SE.mat file
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
%       pointsOn:   Boolean. Determines the plotting type. If pointsOn is
%                   false, use rasters. If pointsOn is true use point
%                   rasters. Point rasters needs some work as MATLAB has
%                   difficulty controlling the size of marker points below
%                   a certain size.
%
%	Outputs:
%		There are no function outputs. But there is a saved file output.
%		Figures
%
%	To Use:
%		Run function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = Rastergram( SpikeData, saveData, channel, stimOn, stimTimes, stimColor, timeBounds, events, pointsOn )

if nargin < 9
    fprintf('Not enough input arguments to Rastergram. Quitting program. DEBUG_nargin\n');
    return
end %END IF

spikeTimes = SpikeData.timeStamps; % Cell containing the spike times for all channels
trialGap   = 1; % "Distance" between trials
H = figure(1);   % Initializes figure

% Ensures that there is at least one arbitrary event
if ~isempty(events)
    numEvents = length(events);
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

for j = channel
    for i = 1:numEvents
        
        startTime = events(i) + timeBounds(1);	% Seconds
        endTime   = events(i) + timeBounds(2);	% Seconds
        
        startTime = startTime * 10e5;			% micro seconds
        endTime   = endTime   * 10e5;			% micro seconds
        
        % Finds the indices of spikes that occur between [startTime] and [endTime]
        trialSpikesIdx = spikeTimes{j} > startTime & spikeTimes{j} < endTime;
        
        % Selects the timeStamps using the indices found above
        trialSpikes{i,1} = spikeTimes{j}(trialSpikesIdx);
        
        % Normalizes the time to 0 seconds
        trialSpikes{i,1} = trialSpikes{i,1} - events(i)*10e5;
        
    end % END FOR
end % END FOR

% Plotting
if ~pointsOn
    hold on
    
    if stimOn
        for q = 1:length(stimTimes/2)
            bgx = [stimTimes(1)-events(1), stimTimes(2)-events(1)];
            bgxx = [bgx, fliplr(bgx)];
        end % END FOR
    end % END IF
    
    trialNum = length(trialSpikes);
    
    bgData = patch(bgxx, [trialNum*11, trialNum*11, 0, 0], 1);
    
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
            
            yy(2:3:3*spikeNum) = (trialNum - (i-1)) * 10 - trialGap;
            yy(1:3:3*spikeNum) = (trialNum - (i-1)-1) * 10;
            
            plot(xx, yy, 'k')
            
        end % END IF
    end % END FOR
    hold off
    
    set(bgData, 'FaceColor', stimColor(1))
    set(bgData, 'EdgeColor', 'none')
    set(bgData, 'FaceAlpha', 0.25)
    
    xlim([xx(1), xx(end)])
    ylim([0, trialNum*10])
    
    set(gca, 'XTick', [0])
    set(gca, 'XTickLabel',{''})
    
    set(gca, 'YTick', [0])
    set(gca, 'YTickLabel',{''})
    
    set(gca, 'Visible', 'Off')
    
    set(gcf, 'Units', 'inches')
    set(gcf, 'Position',[0 0 3 1.5])
    set(gcf, 'PaperUnits','inches','PaperPosition',[0 0 3 1.5])
    
    shortName = ['Raster_', SpikeData.expDate, '_chan', num2str(j)];
    
    if saveData
        print('-dpng', ['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.png'], '-r100');
        saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.fig']);
        %			saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.png']);
        saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.eps'],'epsc2');
    end % END IF
    
elseif pointsOn
        hold on
        
        if stimOn
            for q = 1:length(stimTimes/2)
                bgx = [stimTimes(1)-events(1), stimTimes(2)-events(1)];
                bgxx = [bgx, fliplr(bgx)];
            end % END FOR
        end % END IF
        
        trialNum = length(trialSpikes);
        
        bgData = patch(bgxx, [trialNum*11, trialNum*11, 0, 0], 1);     
        
        for i=1:length(trialSpikes)
            % Checks if there were spikes within the trial region
            if ~isempty(trialSpikes{i})
                spikeNum = length(trialSpikes{i});
                
                xx = ones(1, spikeNum) * nan;
                yy = ones(1, spikeNum) * nan;
                
                % Convert to seconds and organize xx in to sets of 3 points
                xx = double(trialSpikes{i}).*(1/10e5);
                
                yy(1:end) = (i);
                
%                 plot(H, xx, yy, 'k')
                scatter(xx,yy,30,'o','MarkerEdgeColor','none', 'MarkerFaceColor', 'b', 'linewidth', 4)
            end % END IF
        end % END FOR
        hold off
        
        set(bgData, 'FaceColor', stimColor(1))
        set(bgData, 'EdgeColor', 'none')
        set(bgData, 'FaceAlpha', 0.25)
        
        xlim([xx(1), xx(end)])
        ylim([0, i+1])
        
        set(gca, 'XTick', [0])
        set(gca, 'XTickLabel',{''})
        
        set(gca, 'YTick', [0])
        set(gca, 'YTickLabel',{''})
        
        set(gca, 'Visible', 'Off')
        
        set(gcf, 'Units', 'inches')
        set(gcf, 'Position',[1 1 2 1])
        set(gcf, 'PaperUnits','inches','PaperPosition',[0 0 2 1])
        
        shortName = ['Raster_', SpikeData.expDate, '_chan', num2str(j), '_points'];
        
        if saveData
            print('-dpng', ['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.png'], '-r100');
            saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.fig']);
            %			saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.png']);
            saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.eps'],'epsc2');
        end % END IF
        
    else
        fprintf('Something went horribly wrong in Rastergram. DEBUG_Plot\n');
    end % END IF
end % END FUNCTION

% EOF