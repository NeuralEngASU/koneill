%% Load Data

szFilePath = 'D:\PLI\SeizureDetection\ProcessedPLI';

fileNameRaw     = '2014PP04Sz4_PLI_winSize1.mat';
fileNameNotched = '2014PP04Sz4_Filt_PLI_winSize1.mat';
fileName0_50    = '2014PP04Sz4_0_50_PLI_winSize1.mat';
fileName50_100  = '2014PP04Sz4_50_100_PLI_winSize1.mat';
fileName100_150 = '2014PP04Sz4_100_150_PLI_winSize1.mat';
fileName150_200 = '2014PP04Sz4_150_200_PLI_winSize1.mat';
fileName200_250 = '2014PP04Sz4_200_250_PLI_winSize1.mat';

szDataFile = '2014PP04Sz4.mat';

% Load Raw PLI
load(fullfile(['D:\PLI\SeizureDetection\Sz\HilbertFirst'], fileNameRaw));
szRawChanPairNums = chanPairNums;
szRawHeader = Header;
szRawPLI = p;
szRawR = r;
szRawParams = params;
szRawPhi = phi;

% Load Notched PLI
load(fullfile(szFilePath, fileNameNotched));
szNotchChanPairNums = chanPairNums;
szNotchHeader = Header;
szNotchPLI = p;
szNotchR = r;
szNotchParams = params;
szNotchPhi = phi;

% Load 0:50 PLI
load(fullfile(szFilePath, fileName0_50));
sz0ChanPairNums = chanPairNums;
sz0Header = Header;
sz0PLI = p;
sz0R = r;
sz0Params = params;
sz0Phi = phi;

% Load 50:100 PLI
load(fullfile(szFilePath, fileName50_100));
sz50ChanPairNums = chanPairNums;
sz50Header = Header;
sz50PLI = p;
sz50R = r;
sz50Params = params;
sz50Phi = phi;

% Load 100:150 PLI
load(fullfile(szFilePath, fileName100_150));
sz100ChanPairNums = chanPairNums;
sz100Header = Header;
sz100PLI = p;
sz100R = r;
sz100Params = params;
sz100Phi = phi;

% Load 150:200 PLI
load(fullfile(szFilePath, fileName150_200));
sz150ChanPairNums = chanPairNums;
sz150Header = Header;
sz150PLI = p;
sz150R = r;
sz150Params = params;
sz150Phi = phi;

% Load 200:250 PLI
load(fullfile(szFilePath, fileName200_250));
sz200ChanPairNums = chanPairNums;
sz200Header = Header;
sz200PLI = p;
sz200R = r;
sz200Params = params;
sz200Phi = phi;

% Clear old data
clear chanPairNums Header p params phi r

% Load Raw Data
load(['E:\data\human CNS\EMD\Sz\clips\', szDataFile])
data = detrend(data');

%%

refChan = 65; % 2014PP04Sz4
% refChan = 65; % 2012PP05Sz1

badChan = []; % For 2014PP04Sz4
% badChan = [33,32]; % For 2012PP05Sz1
chanLoop = 1:64;

for ii = 1:length(badChan)
    chanLoop = chanLoop(chanLoop ~= badChan(ii));
end % END FOR

winSize = 5;
winSizeAdj = winSize-1;

timeWindowIdxBounds = [[1:600-winSizeAdj];[1+winSizeAdj:600]];

integralDiff = [];

% timeWindowBounds = [[-5:0.25:4.75]; [-4.75:0.25:5]];

szDataNames = {'szRawPLI','szNotchPLI','sz0PLI','sz50PLI','sz100PLI', 'sz150PLI', 'sz200PLI'};
szPairNames = {'szRawChanPairNums','szNotchChanPairNums','sz0ChanPairNums','sz50ChanPairNums','sz100ChanPairNums', 'sz150ChanPairNums', 'sz200ChanPairNums'};


for type = 1:7
    
    eval(['szPLI = ', szDataNames{type}, ';']);
    eval(['szChanPairNums = ', szPairNames{type}, ';']);
    
    for chanIdx = chanLoop
        
        % Find Channels to Plot (3x3 grid around chanIdx)
        gridDim = [8,8];
        [ii, jj] = ind2sub(gridDim, chanIdx);
        
        iiNew = [ii-1, ii, ii+1];
        jjNew = [jj-1, jj, jj+1];
        
        % Find index pairs
        [idx2, idx1] = find(true(numel(iiNew),numel(jjNew)));
        indPairs = [reshape(iiNew(idx1), [], 1), reshape(jjNew(idx2), [], 1)];
        
        % Remove invalid pairs (where there is a 0 index)
        indParisIdx = (indPairs(:,1) > 0) & (indPairs(:,2) > 0) & (indPairs(:,1) <= gridDim(1)) & (indPairs(:,2) <= gridDim(2));
        indPairs2 = indPairs(indParisIdx,:);
        
        % Convert index notation to channel number
        chansPlot = sort(sub2ind([8,8], indPairs2(:,1), indPairs2(:,2)), 'ascend');        
        desiredRef = [chanIdx];


%         desiredRef = refChan;
%         chansPlot = chanIdx;
        
        for ii = 1:length(badChan)
            chansPlot = chansPlot(chansPlot ~= badChan(ii));
        end % END FOR
        
        % Define channel pairs
        if isempty(desiredRef)
            desiredChanPairs = nchoosek(sort(unique(chansPlot),'ascend'),2);
        else
            desiredChanPairs = [repmat(desiredRef, size(chansPlot(:),1), 1), chansPlot(:)];
        end % END IF isempty(desiredRef)
        
        for ii = 1:size(desiredChanPairs,1)
            if desiredChanPairs(ii,1) > desiredChanPairs(ii,2)
                desiredChanPairs(ii,:) = fliplr(desiredChanPairs(ii,:));
            end % END IF
        end
        
        idx = zeros(size(szChanPairNums,1),1);
        % Find idicies
        for jj = 1:size(desiredChanPairs,1)
            idx = idx | ((szChanPairNums(:,1) == desiredChanPairs(jj,1)) & (szChanPairNums(:,2) == desiredChanPairs(jj,2)));
        end
        
        idx = find(idx==1);
        
%         sizePLI = size(szPLI,1) * ( size(szPLI,1) <  size(nonSzPLI,1)) + size(nonSzPLI,1) * ( size(szPLI,1) >=  size(nonSzPLI,1));
        sizePLI = size(szPLI,1);
        
%         varSzPLI = var(szPLI(:,idx)');
%         varNonSzPLI = var(nonSzPLI(:,idx)');
        
        if length(idx)==1
            meanSzPLI = szPLI(1:sizePLI,idx)';
%             meanNonSzPLI = nonSzPLI(1:sizePLI,idx)';
        else
            meanSzPLI = mean(szPLI(1:sizePLI,idx)');
%             meanNonSzPLI = mean(nonSzPLI(1:sizePLI,idx)');
        end
        
        timeMax = size(szPLI,1);
        sizeMax = 500*timeMax;
        
        time = linspace(-5,5,sizePLI);
        
        % Channel, time window, sz/nonSz
%         rawIntegral(chanIdx,1,1) = sum(meanSzPLI)/sizePLI;
%         rawIntegral(chanIdx,1,2) = sum(meanNonSzPLI)/sizePLI;
        
        for tw = 1:size(timeWindowIdxBounds,2)
            
            timeWindowIdxBoundsMin = timeWindowIdxBounds./60-5;
            
            %     timeWindow = time >= timeWindowIdxBoundsMin(1,tw) & time < timeWindowIdxBounds(2,tw);
            timeWindow = timeWindowIdxBounds(1,tw) : timeWindowIdxBounds(2,tw);
            %     timeSize = 60* diff(timeWindowIdxBoundsMin(:,tw));
            timeSize = length(timeWindow);
            
            
            % Channel, time window
            %     integralDiff(chanIdx,tw) = abs(sum(meanSzPLI(timeWindow))/timeSize - sum(meanNonSzPLI(timeWindow))/timeSize);
            integralDiff(chanIdx,tw,type) = abs(sum(meanSzPLI(timeWindow))/timeSize);
            
        end % END FOR each timeWindow
    end % END FOR
end % END FOR type

%% Plot Test
figure;

testChan = 1;
subplot(7,1,1)
plot(integralDiff(testChan,:,1))
ylim([0,1])

subplot(7,1,2)
plot(integralDiff(testChan,:,2))
ylim([0,1])

subplot(7,1,3)
plot(integralDiff(testChan,:,3))
ylim([0,1])

subplot(7,1,4)
plot(integralDiff(testChan,:,4))
ylim([0,1])

subplot(7,1,5)
plot(integralDiff(testChan,:,5))
ylim([0,1])

subplot(7,1,6)
plot(integralDiff(testChan,:,6))
ylim([0,1])

subplot(7,1,7)
plot(integralDiff(testChan,:,7))
ylim([0,1])

%%

figure;
[Pxy,F] = mscohere(data(:,1),data(:,2),hamming(100),80,100,Fs);

plot(F,Pxy)
title('Magnitude-Squared Coherence')
xlabel('Frequency (Hz)')
grid

% EOF