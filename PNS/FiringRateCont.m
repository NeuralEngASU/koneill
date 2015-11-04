function [ rateMat ] = FiringRateCont( trialStruct, trialType, trials, channels, trialLen, spikeTimes )

numIn = nargin;

%%%%% INIT %%%%%
fs         = 30000; % default sampling rate
flagLabel  = 1;     % boolean to show labels
spikeWidth = 1;     % spike thickness
spikeColor = 'k';   % spike color
trialGap   = 1.5;   % distance between trials
buffer     = 15000; % beginning and end buffer in samples

%%%%% MAIN %%%%%
trialSize = trialLen + 2* buffer;
rasterTrials = eval(['trialStruct.', trialType,'(trials,:)']);
numTrials = length(trials);
numChans     = size(channels, 1);
rateMat = zeros(numChans, trialSize, numTrials);
boxWin = BoxcarWindow(300, 30000);

% loops over channels
for j = 1:numChans
    % Loops over the number of trials
    for k = 1:numTrials
        
        % Size of each trial to plot (Samples)
        trialSize = trialLen + 2* buffer;
        startTime = rasterTrials(k,4) - buffer;
%         timeMat = startTime:(startTime + trialSize);
        
        trialSpikesIdx = spikeTimes{channels(j)} > startTime & spikeTimes{channels(j)} < [startTime + trialSize];
        trialSpikes = spikeTimes{channels(j)}(trialSpikesIdx);
        
        numSpike = length(trialSpikes);
        if numSpike
            % Makes spikes impulses
            rateMat(j, [trialSpikes - startTime + 1], k) = 1;
            
            % Convolves impulses with window
            tmpMat = conv(rateMat(j, :, k), boxWin);
            rateMar(j, :, k) = tmpMat(1:end-9000);
            
            % low pass filter on convolution
            T   = 1/fs; % sec
            tau = 0.5;  % sec
            a = T/tau;
            rateMat(j, :, k) = filter(a, [1, a-1], rateMat(j, :, k));
        end
    end % END FOR
end % END FOR



end % END FUNCTION

