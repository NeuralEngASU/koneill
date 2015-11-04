function [L] = ClassifyPCA( trialStruct, spikeTimes, PCA, sampleTime, currTime )


numIn = nargin;

%%%%% INIT %%%%%
dataFs = 30000; % default sampling rate

%%%%% MAIN %%%%%
numChans     = size(PCA.Electrodes, 1);
currTime     = currTime; % Current time in the Data, samples
PCA.Channels = e2c(PCA.Electrodes,'pns');


[V,D] = eig(cov(PCA.Results.TrainMat));
V = fliplr(V); D = fliplr(D); %organizing vectors so most variance is in left column
Var90 = find(cumsum(sum(D))/sum(D(:))<=0.9,1,'last'); %finding first set of vectors with 90 percent of the variance
if isempty(Var90)
    Var90 = 1;
end

boxwin = 0.3; %boxcar window in sec
RateDS.WfTS = cell(96,1);
RateDS.Data = zeros(96,15*sampleTime); %total length for 15 S/sec
Hd = design(fdesign.lowpass('N,F3dB',4,7.5,600),'butter'); %spikes initially downsampled to 600S/sec for convolution, then down to 15
for k = 1:numChans
%     disp(k)
    
    trialSpikesIdx = spikeTimes{PCA.Channels(k)} > currTime - 2*dataFs & spikeTimes{PCA.Channels(k)} < currTime;
    
    trialSpikes = spikeTimes{PCA.Channels(k)}(trialSpikesIdx);
    
    TS = double(trialSpikes);
    RateDS.WfTS(k) = {TS./2000}; % timestamps downsampled to 15 S/sec
    RateB = zeros(1,15*sampleTime*40); % total length for 600 S/sec
    if ~isempty(TS)
        RateB(ceil((TS - (currTime-2*dataFs))./50)) = 1; %timestamps at 600 S/sec
        RateB = conv2(RateB,ones(1,boxwin*600)./boxwin);
        RateB = RateB(1:15*sampleTime*40);
        RateB = filter(Hd,RateB); %filter at 7.5Hz
    end % END IF
    RateDS.Data(k,:) = RateB(1:40:end); %downsample from 600 S/sec to 15 S/sec
end % END FOR
clear 'RateB'

testMat = RateDS.Data(PCA.Channels, :);

testMat = reshape(testMat, 1, []);

PCTestMat = testMat*V(:,1:Var90);

L = classify(PCTestMat,PCA.Results.PCTrainMat,PCA.Results.TrainMatID(:,1)); %classifying


end % END FUNCTION

