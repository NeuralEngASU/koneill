function [ rateMat ] = FiringRatePSTH( trialStruct, trialType, trials, channels, trialLen, spikeTimes  )


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
% rateMat = zeros(numChans, trialSize, numTrials);
boxWin = BoxcarWindow(300, 30000);

% loops over channels
for j = 1:numChans
    % Loops over the number of trials
    for k = 1:numTrials
        
        % Size of each trial to plot (Samples)
%         trialSize = trialLen + 2* buffer;
        trialSize = trialLen;
        startTime = rasterTrials(k,4);
%         timeMat = startTime:(startTime + trialSize);
        
        trialSpikesIdx = spikeTimes{channels(j)} > startTime - 2*fs & spikeTimes{channels(j)} < [startTime + trialSize + 2*fs];
        trialSpikes = spikeTimes{channels(j)}(trialSpikesIdx);
        
        trialSpikes = trialSpikes - startTime;
        
        AP.spiketimes = double(trialSpikes)./fs;
%         eval(['AP.Trial', num2str(k), ' = double(trialSpikes)./fs;']);
        
        
        numSpike = length(trialSpikes);
        if numSpike
            
%             h = figure(1);
%             hold on
            [R, t, E] = psth(AP, 0.300, 'r', [0, 3], 1);
%             close(h)
%             hold off
            
            sz = size(R);
            
            x=linspace(0, 3, length(R));

            if sz(1) < sz(2)
                R = R';
                t = t';
                E = E';
            end;
            
            RCell{k} = R;
            tCell{k} = t;
            ECell{k} = E;
%             
%             hold on
%             figure (2)
%             plot(x, R, 'b')
%             hold off
        else
            RCell{k} = NaN;
            tCell{k} = NaN;
            ECell{k} = NaN;
        end
        
%         
%         numSpike = length(trialSpikes);
%         if numSpike
%             % Makes spikes impulses
%             rateMat(j, [trialSpikes - startTime + 1], k) = 1;
%             
%             % Convolves impulses with window
%             tmpMat = conv(rateMat(j, :, k), boxWin);
%             rateMar(j, :, k) = tmpMat(1:end-9000);
%             
%             % low pass filter on convolution
%             T   = 1/fs; % sec
%             tau = 0.5;  % sec
%             a = T/tau;
%             rateMat(j, :, k) = filter(a, [1, a-1], rateMat(j, :, k));
%         end
    end % END FOR
    
    count = 1;
    for c = 1:length(RCell)
        
        if ~isnan(RCell{c})
            RMat(:,count) = RCell{c};
            EMat(:,count) = ECell{c};
            count = count+1;
        end
        
    end % END FOR
    
    
    
    RMat = mean(RMat, 2);
    EMat = mean(EMat, 2);
    
    if j == 1
        rateMat = zeros(length(RMat), numChans, 2);
    end;
    rateMat(:, j, 1) = RMat;
    rateMat(:, j, 2) = EMat;
    
    
end % END FOR


end % END FUNCTION

% EOF