function [ firingRate ] = FiringRate( trialStruct, trialType, trials, channels, trialLen, spikeTimes )

numIn = nargin;

%%%%% INIT %%%%%
fs         = 30000; % default sampling rate
flagLabel  = 1;     % boolean to show labels
spikeWidth = 1;     % spike thickness
spikeColor = 'k';   % spike color
trialGap   = 1.5;   % distance between trials
buffer     = 15000; % beginning and end buffer in samples

%%%%% MAIN %%%%%

rasterTrials = eval(['trialStruct.', trialType,'(trials,:)']);
numTrials = length(trials);
rasterMat = zeros(1, numTrials);
numChans     = size(channels, 1);

rasterMat = zeros(numChans, numTrials);

% loops over channels
for j = 1:numChans
    % Loops over the number of trials
    for k = 1:numTrials
        
        % Size of each trial to plot (Samples)
        trialSize = trialLen;
        startTime = rasterTrials(k,4);
        
        trialSpikesIdx = spikeTimes{channels(j)} > startTime & spikeTimes{channels(j)} < [startTime + trialSize];
        trialSpikes = spikeTimes{channels(j)}(trialSpikesIdx);
        
        spikeNum = length(trialSpikes);
        
        rasterMat(j,k) = spikeNum;
    end % END FOR
end % END FOR

firingRate = rasterMat./(trialLen/fs);

firingRate = firingRate'; % Transposed for boxplot() function.

end % END FUNCTION

