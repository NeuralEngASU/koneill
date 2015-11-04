function [] = RasterPlotCorr( trialStruct, trialType, trials, corrStruct, trialLen, spikeTimes, H, varargin )

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
buffer     = 7500;  % beginning and end buffer in samples

switch numIn
    case 7 % Use defaults on Spike Times
        if isempty(H)
            figure;
            H=gca;
        end
    case 8 % Alters sampling rate
        figure;
        fs = varargin{1};
    case 9 % Alters labels
        figure;
        fs        = varargin{1};
        flagLabel = varargin{2};
    case 10 % Spike Width
        figure;
        fs         = varargin{1};
        flagLabel  = varargin{2};
        spikeWidth = varargin{3};
    case 11 % Spike Color
        figure;
        fs         = varargin{1};
        flagLabel  = varargin{2};
        spikeWidth = varargin{3};
        spikeColor = varargin{4};
    otherwise
        error ('Invalid Arguments');
end

%%%%% MAIN %%%%%

corrChans    = eval(['corrStruct.', trialType,'(:,:)']);
rasterTrials = eval(['trialStruct.', trialType,'(trials,:)']);
numTrials    = length(trials);
numChans     = size(corrChans, 1);

hold on
% Loops over the number of trials
for i = 1:numChans
    for k = 1:numTrials
        
        % Size of each trial to plot (Samples)
        trialSize = trialLen + 2* buffer;
        startTime = rasterTrials(k,4) - buffer;
        
        trialSpikesIdx = spikeTimes{corrChans(i,2)} > startTime & spikeTimes{corrChans(i,2)} < [startTime + trialSize];
        trialSpikes = spikeTimes{corrChans(i,2)}(trialSpikesIdx);
        
        if ~isempty(trialSpikes)
            spikeNum = length(trialSpikes);
            
            xx = ones(1, 3 * spikeNum) * nan;
            yy = ones(1, 3 * spikeNum) * nan;
            
            % Convert to seconds and organize xx in to sets of 3 points
            xx(2:3:3*spikeNum) = double((trialSpikes - startTime)).*(1/fs);
            xx(1:3:3*spikeNum) = xx(2:3:3*spikeNum) - eps(1);
            xx(3:3:3*spikeNum) = xx(2:3:3*spikeNum) + eps(1);
            
            yy(2:3:3*spikeNum) = (numChans - (i-1)) * 10 - trialGap;
            yy(1:3:3*spikeNum) = (numChans - (i-1)-1) * 10;
            
            plot(H, (xx - (buffer/fs) + (k-1)*trialSize/fs), yy, spikeColor, 'linewidth', spikeWidth)
           
        end % END IF
    end % END FOR
end % END FOR

xFlag = ones(1, 6*numTrials) * nan;
yFlag = ones(1, 6*numTrials) * nan;

xFlag(1:6:6*numTrials) = ([1:numTrials]-1) .* (trialSize) ./fs;
xFlag(4:6:6*numTrials) = ([1:numTrials]) .* (trialSize) ./fs - 2*buffer/fs;
xFlag(2:6:6*numTrials) = xFlag(1:6:6*numTrials) + eps(1);
xFlag(5:6:6*numTrials) = xFlag(4:6:6*numTrials) + eps(1);

yFlag(1:6:6*numTrials) = 0;
yFlag(4:6:6*numTrials) = 0;
yFlag(2:6:6*numTrials) = numChans*10;
yFlag(5:6:6*numTrials) = numChans*10;

plot(H, xFlag, yFlag, 'b:')

hold off

xlim([-buffer/fs, (numTrials*trialSize)/fs])
ylim([0, numChans*10])

if flagLabel
    ylabel('Channels')
    xlabel(sprintf('Trials (%d sec/trial)', trialLen/fs))
end % END IF

set(gca, 'ytick', [1:numChans] * 10 - 5)
set(gca, 'yticklabel', [flipud(corrChans(:,2))])
set(gca, 'xtick', [(trialSize/2 - buffer)  + trialSize*([1:numTrials]-1)]./fs)
set(gca, 'xticklabel', [1:numTrials])

end % END FUNTION

% EOF
