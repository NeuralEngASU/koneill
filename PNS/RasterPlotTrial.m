function [rasterMat] = RasterPlotTrial( trialStruct, trialType, trials, trialLen, spikeTimes, H, varargin )

% trialStruct
% Trial Type
% Trials
% Spike Times
% Varargin

numIn = nargin;

%%%%% INIT %%%%%
fs         = 30000; % default sampling rate
flagLabel  = 1;     % boolean to show labels
spikeWidth = 1;     % spike thickness
spikeColor = 'k';   % spike color
trialGap   = 1.5;   % distance between trials
buffer     = 15000; % beginning and end buffer in samples

switch numIn
    case 6 % Use defaults on Spike Times
        if isempty(H)
            figure;
            H=gca;
        end
    case 7 % Alters sampling rate
        figure;
        fs = varargin{1};
    case 8 % Alters labels
        figure;
        fs        = varargin{1};
        flagLabel = varargin{2};
    case 9 % Spike Width
        figure;
        fs         = varargin{1};
        flagLabel  = varargin{2};
        spikeWidth = varargin{3};
    case 10 % Spike Color
        figure;
        fs         = varargin{1};
        flagLabel  = varargin{2};
        spikeWidth = varargin{3};
        spikeColor = varargin{4}; 
    otherwise
        error ('Invalid Arguments');
end

%%%%% MAIN %%%%%

rasterTrials = eval(['trialStruct.', trialType,'(trials,:)']);
numTrials = length(trials);
rasterMat = zeros(1, numTrials);

hold on
% Loops over the number of trials
for k = 1:numTrials
    
    % Size of each trial to plot (Samples)
    trialSize = trialLen + 2* buffer;
    startTime = rasterTrials(k,4) - buffer;
    
    trialSpikesIdx = spikeTimes > startTime & spikeTimes < [startTime + trialSize];
    trialSpikes = spikeTimes(trialSpikesIdx);
    
    if ~isempty(trialSpikes)
        spikeNum = length(trialSpikes);
        
        xx = ones(1, 3 * spikeNum) * nan;
        yy = ones(1, 3 * spikeNum) * nan;
        
        % Convert to ms and organize xx in to sets of 3 points
        xx(2:3:3*spikeNum) = (trialSpikes - startTime).*(1000/fs);
        xx(1:3:3*spikeNum) = xx(2:3:3*spikeNum) - eps(1);
        xx(3:3:3*spikeNum) = xx(2:3:3*spikeNum) + eps(1);
        
        yy(2:3:3*spikeNum) = (numTrials - (k-1)) * 10 - trialGap;
        yy(1:3:3*spikeNum) = (numTrials - (k-1)-1) * 10;
        
        plot(H, xx - (buffer*1000/fs), yy, spikeColor, 'linewidth', spikeWidth)
        
        rasterMat(k) = spikeNum;
    end % END IF
        
end % END FOR

hold off

xlim([-buffer*1000/fs, (trialSize)*1000/fs])
ylim([0, numTrials*10])

if flagLabel
    ylabel('Trial')
    xlabel('Time, ms')
end % END IF

set(gca, 'ytick', [])

end % END FUNTION

% EOF
