function [  ] = InterSpikeInterval( input_args )

% Defaults
CHANNELS = 1;
NUMBINS = 100;
MAXTIME = 0.5;
BINCOLOR = 'b';
NAMEMOD = '';

% Parse varargin
for i = 1:2:length(varargin)
    eval([upper(varargin{i}) '=varargin{i+1};']);
end

% Rename Variables to be consistant with convention
numBins = NUMBINS;
maxTime = MAXTIME;
binColor = BINCOLOR;
nameMod = NAMEMOD;
channels = CHANNELS;

% Allocate channelIdx
channelIdx = zeros(length(SpikeChans));

for i = 1:length(channels);

    % Find all spikes from selected channel
    channelIdx = SpikeChan == channels(i);
    
    % Calculate the difference in time between the spikes
    timeDiff = diff(SpikeTs(channelIdx));
    
    % Sort intervals in ascending order.
    timeDiff = sort(timeDiff, 'ascend');
    
    % Remove all intervals greater than maxTime
    timeDiff = timeDiff(timeDiff<=maxTime);
    
    % Calculate bin placement
    binLoc = maxTime/numBins;
    binCenter = binLoc/2;
    binCenters = linspace(binCenter, maxTime-binCenter, numBins);
    
    % Plot histogram
    hist(timeDiff, binCenters)
    
%     if ischar(binColor)
%         binColor = ColorFinder(binColor);
%     end % END IF
    
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor', [0,0,1])
    %set(h,'EdgeColor', 'none')
    
end % END FOR
end % END FUNCTION

% EOF