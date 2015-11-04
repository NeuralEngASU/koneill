function [ tuneC ] = TuningCurve( trialStruct, trials, channels, trialLen, spikeTimes  )

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
numTrials = length(trials);
numChans     = size(channels, 1);
% rateMat = zeros(numChans, trialSize, numTrials);
boxWin = BoxcarWindow(300, 30000);

% fieldList = fieldnames(trialStruct);

fieldList = {'T20000000', 'T02000000', 'T00200000', 'T00020000', 'T00002000', ... % TIMRL Flex
             'T00000200', 'T00000020', 'T00000002', ...                           % IRL Abduction.
             'T10000000', 'T01000000', 'T00100000', 'T00010000', 'T00001000', ... % TIMRL Extenstion
             'T00000000'};                                                        % Null state

tuneC = zeros(numChans, numTrials);

% Loops over movement types
for i = 1:length(fieldList)
    
    if i ~= 14
        rasterTrials = eval(['trialStruct.', fieldList{i},'(trials,:)']);
    else
        rasterTrials = eval(['trialStruct.', fieldList{1},'(trials,:)']);
    end
    % loops over channels
    for j = 1:numChans
        % Loops over the number of trials
        for k = 1:numTrials
            
                        
            % Size of each trial to plot (Samples)
            %         trialSize = trialLen + 2* buffer;
            trialSize = trialLen;
            startTime = rasterTrials(k,4);
            %         timeMat = startTime:(startTime + trialSize);
            
            % If movement is not null state, average between 0.5 and 2.5 
            % sec of trial
            if i ~= 14
                trialSpikesIdx = spikeTimes{channels(j)} > startTime + fs/2 & spikeTimes{channels(j)} < [startTime + fs/2 + 2*fs];
            else % Average over the previous 2 seconds is movement is null
                trialSpikesIdx = spikeTimes{channels(j)} > startTime - 2*fs & spikeTimes{channels(j)} < startTime;
            end
            trialSpikes = spikeTimes{channels(j)}(trialSpikesIdx);
            
            trialSpikes = trialSpikes - startTime;
            
            tuneC(j,k,i) = length(trialSpikes)/2;
            
        end % END FOR
    end % END FOR
end % END FOR
end % END FUNTION

% EOF

