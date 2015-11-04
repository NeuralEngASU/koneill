%% Animated Integral Plots of 3x3 Grid 

szFilePath = 'D:\PLI\SeizureDetection\Sz\HilbertFirst';
nonSzFilePath = 'D:\PLI\SeizureDetection\NonSz\HilbertFirst';

szFileName = dir([szFilePath, '\*.mat']);
nonSzFileName = dir([nonSzFilePath, '\*.mat']);
numFiles = size(szFileName,1);

fileOfInterest = 1;

% szFileName(fileOfInterest).name = '2012PP05Sz1_PLI_winSize1.mat';
% nonSzFileName(fileOfInterest).name = '2012PP05NonSz1_PLI_winSize1.mat';
% szDataFile = '2012PP05Sz1.mat';

szFileName(fileOfInterest).name = '2014PP04Sz4_PLI_winSize1.mat';
nonSzFileName(fileOfInterest).name = '2014PP04NonSz4_PLI_winSize1.mat';
szDataFile = '2014PP04Sz4.mat';

% Load Sz
% load(fullfile(szFilePath, szFileName(fileOfInterest).name));
load(fullfile(szFilePath, szFileName(fileOfInterest).name));
szChanPairNums = chanPairNums;
szHeader = Header;
szPLI = p;
szR = r;
szParams = params;
szPhi = phi;

% Load NonSz
% load(fullfile(nonSzFilePath, nonSzFileName(fileOfInterest).name));
load(fullfile(nonSzFilePath, nonSzFileName(fileOfInterest).name));

nonSzChanPairNums = chanPairNums;
nonSzHeader = Header;
nonSzPLI = p;
nonSzR = r;
nonSzParams = params;
nonSzPhi = phi;
% Clear old data
clear chanPairNums Header p params phi r

integralDiff = [];

timeWindowBounds = [[-5:0.25:4.75]; [-4.75:0.25:5]];

for chanIdx = 1:64

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

sizePLI = size(szPLI,1) * ( size(szPLI,1) <  size(nonSzPLI,1)) + size(nonSzPLI,1) * ( size(szPLI,1) >=  size(nonSzPLI,1));

varSzPLI = var(szPLI(:,idx)');
varNonSzPLI = var(nonSzPLI(:,idx)');

meanSzPLI = mean(szPLI(1:sizePLI,idx)');
meanNonSzPLI = mean(nonSzPLI(1:sizePLI,idx)');

timeMax = size(szPLI,1);
sizeMax = szHeader.params.Fs*timeMax;

time = linspace(-5,5,sizePLI);

% Channel, time window, sz/nonSz
rawIntegral(chanIdx,1,1) = sum(meanSzPLI)/sizePLI;
rawIntegral(chanIdx,1,2) = sum(meanNonSzPLI)/sizePLI;

for tw = 1:size(timeWindowBounds,2)
    
    timeWindow = time >= timeWindowBounds(1,tw) & time < timeWindowBounds(2,tw);
    timeSize = 60* diff(timeWindowBounds(:,tw));
    
    % Channel, time window
    integralDiff(chanIdx,tw) = abs(sum(meanSzPLI(timeWindow))/timeSize);
    
end % END FOR each timeWindow
end % END FOR

%% Animate

%%

% integralDiff = integralDiff./max(integralDiff(:));

load(['E:\data\human CNS\EMD\Sz\clips\', szDataFile])

figure;
mapCol = gray(128);
border = 0;

layout = reshape(1:64,8,8);


for kk = 1:size(integralDiff,2)
    clf
    ax1 = subplot(1,2,1);
    hold on
    
    plotData = detrend(data(1,1:sizePLI*500));
    
    patchx = [timeWindowBounds(1,kk), timeWindowBounds(2,kk), timeWindowBounds(2,kk), timeWindowBounds(1,kk)];
    patchy = [2*min(plotData), 2*min(plotData), 2*max(plotData), 2*max(plotData)];
    
    pData = patch(patchx, patchy,1);
    
    set(pData, 'EdgeColor', 'none')
    set(pData, 'FaceColor', [0.5,0.5,0.5])
    set(pData, 'FaceAlpha', 0.25)
    
    
    
    plot(linspace(-5,5,sizePLI*500), plotData)
    title('Detrended Data from Channel 1')
    xlabel('Time, minutes')
    ylabel('Voltage, uV')
    ylim([1.1*min(plotData), 1.1*max(plotData)])
    
    ax2 = subplot(1,2,2);
    % Find Lower Left corner
    chans = 1:64;
    [~,chanLLIdx] = intersect(layout,chans);
    
    % Subdivide each section into [x,y] sections. Where x and y are layout
    % dimensions.
    rowsFix = linspace(0+border,1-border,size(layout,1)+1);
    
    colsFix = linspace(0+border,1-border,size(layout,2)+1);
    
    for jj = 1:length(chanLLIdx)
        
        chanIdx = layout(chanLLIdx(jj));
        [offx, offy] = ind2sub([size(layout,1),size(layout,2)],find(layout==chanIdx));
        
        xSubPos = [rowsFix(1,1), rowsFix(1,end), rowsFix(1,end), rowsFix(1,1)  ] + offx;
        ySubPos = [colsFix(1,1), colsFix(1,1)  , colsFix(1,end), colsFix(1,end)] + offy;
        
        chanIdx2 = layout(chanLLIdx(jj));
        
        [xpos, ypos] = ind2sub([size(layout,1),size(layout,2)],find(layout==chanIdx2));
        
        colorPatch = mapCol(floor(integralDiff(jj,kk) * 127+1),:);
        
        pData  = patch( xSubPos, ySubPos, colorPatch);
        
        
        xlim([0.5,size(layout,1)+1.5])
        ylim([0.5,size(layout,2)+1.5])
        
        set(pData, 'EdgeColor', 'none')
        axis square
        
        titleString = {'Pre-Ictal', 'Ictal', 'Early Post-Ictal', 'Late Post-Ictal'};
        
%         title([titleString{ii}, ': ', num2str(timeWindowBounds(1,ii)), ' to ', num2str(timeWindowBounds(2,ii)), ' minutes'])
        
        set(gca, 'xtick', [1.5:8.5])
        set(gca, 'xticklabel', [1:8])
        
        set(gca, 'ytick', [1.5:8.5])
        set(gca, 'yticklabel', [1:8:57])
        
        xlabel('Channel')
        ylabel('Channel')
        box on
        
    end % END FOR num LL corners
    set(gcf, 'units', 'inches')
    set(gcf, 'Position', [2    2   12    4])
    set(gcf, 'PaperPosition', [2    2   12    4])
    
    pos1=get(ax1, 'Position');
    axes(ax1)
    set(ax1, 'Position', [.055 pos1(2) .575 pos1(4)])
%     
    pos2=get(ax2, 'Position');
    axes(ax2)
    set(ax2, 'Position', [.675 pos2(2) pos2(3) pos2(4)])
    
    drawnow
    
    tmpNumber = num2str(kk);
    if length(tmpNumber) == 1;
        tmpNumber = ['00', tmpNumber];
    elseif length(tmpNumber) == 2
        tmpNumber = ['0', tmpNumber];
    end
    
    print([szDataFile, '_15Sec_Frame_', tmpNumber, '.png'],'-dpng')
    
end

%%

szFilePath = 'D:\PLI\SeizureDetection\Sz\HilbertFirst';
nonSzFilePath = 'D:\PLI\SeizureDetection\NonSz\HilbertFirst';

szFileName = dir([szFilePath, '\*.mat']);
nonSzFileName = dir([nonSzFilePath, '\*.mat']);
numFiles = size(szFileName,1);

fileOfInterest = 1;

% szFileName(fileOfInterest).name = '2012PP05Sz1_PLI_winSize1.mat';
% nonSzFileName(fileOfInterest).name = '2012PP05NonSz1_PLI_winSize1.mat';
% szDataFile = '2012PP05Sz1.mat';

szFileName(fileOfInterest).name = '2014PP04Sz4_PLI_winSize1.mat';
nonSzFileName(fileOfInterest).name = '2014PP04NonSz4_PLI_winSize1.mat';
szDataFile = '2014PP04Sz4.mat';

% Load Sz
% load(fullfile(szFilePath, szFileName(fileOfInterest).name));
load(fullfile(szFilePath, szFileName(fileOfInterest).name));
szChanPairNums = chanPairNums;
szHeader = Header;
szPLI = p;
szR = r;
szParams = params;
szPhi = phi;

% Load NonSz
% load(fullfile(nonSzFilePath, nonSzFileName(fileOfInterest).name));
load(fullfile(nonSzFilePath, nonSzFileName(fileOfInterest).name));

nonSzChanPairNums = chanPairNums;
nonSzHeader = Header;
nonSzPLI = p;
nonSzR = r;
nonSzParams = params;
nonSzPhi = phi;
% Clear old data
clear chanPairNums Header p params phi r

% Load Raw Data
load(['E:\data\human CNS\EMD\Sz\clips\', szDataFile])

%% RefChan

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

timeWindowIdxBounds = [[300:420-winSizeAdj];[300+winSizeAdj:420]];

integralDiff = [];

% timeWindowBounds = [[-5:0.25:4.75]; [-4.75:0.25:5]];

for chanIdx = chanLoop

% Find Channels to Plot (3x3 grid around chanIdx)
% gridDim = [8,8];
% [ii, jj] = ind2sub(gridDim, chanIdx);
% 
% iiNew = [ii-1, ii, ii+1];
% jjNew = [jj-1, jj, jj+1];
% 
% % Find index pairs
% [idx2, idx1] = find(true(numel(iiNew),numel(jjNew))); 
% indPairs = [reshape(iiNew(idx1), [], 1), reshape(jjNew(idx2), [], 1)];
%     
% % Remove invalid pairs (where there is a 0 index)
% indParisIdx = (indPairs(:,1) > 0) & (indPairs(:,2) > 0) & (indPairs(:,1) <= gridDim(1)) & (indPairs(:,2) <= gridDim(2));
% indPairs2 = indPairs(indParisIdx,:);
% 
% % Convert index notation to channel number
% chansPlot = sort(sub2ind([8,8], indPairs2(:,1), indPairs2(:,2)), 'ascend');

chansPlot = chanIdx;

% desiredRef = [chanIdx];
desiredRef = refChan;

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

sizePLI = size(szPLI,1) * ( size(szPLI,1) <  size(nonSzPLI,1)) + size(nonSzPLI,1) * ( size(szPLI,1) >=  size(nonSzPLI,1));

varSzPLI = var(szPLI(:,idx)');
varNonSzPLI = var(nonSzPLI(:,idx)');

if length(idx==1)
    meanSzPLI = szPLI(1:sizePLI,idx)';
    meanNonSzPLI = nonSzPLI(1:sizePLI,idx)';
else
    meanSzPLI = mean(szPLI(1:sizePLI,idx)');
    meanNonSzPLI = mean(nonSzPLI(1:sizePLI,idx)');
end

timeMax = size(szPLI,1);
sizeMax = szHeader.params.Fs*timeMax;

time = linspace(-5,5,sizePLI);

% Channel, time window, sz/nonSz
rawIntegral(chanIdx,1,1) = sum(meanSzPLI)/sizePLI;
rawIntegral(chanIdx,1,2) = sum(meanNonSzPLI)/sizePLI;

for tw = 1:size(timeWindowIdxBounds,2)
    
    timeWindowIdxBoundsMin = timeWindowIdxBounds./60-5;
    
%     timeWindow = time >= timeWindowIdxBoundsMin(1,tw) & time < timeWindowIdxBounds(2,tw);
    timeWindow = timeWindowIdxBounds(1,tw) : timeWindowIdxBounds(2,tw);
%     timeSize = 60* diff(timeWindowIdxBoundsMin(:,tw));
    timeSize = length(timeWindow);
    
    % Channel, time window
%     integralDiff(chanIdx,tw) = abs(sum(meanSzPLI(timeWindow))/timeSize - sum(meanNonSzPLI(timeWindow))/timeSize);
    integralDiff(chanIdx,tw) = abs(sum(meanSzPLI(timeWindow))/timeSize);
    
end % END FOR each timeWindow
end % END FOR

%% Animate and Color (Use SfN_PosterPlots line ~625)

integralDiff2 = integralDiff;
integralDiff = integralDiff./max(integralDiff(:));

for tmpFrame = 1:size(integralDiff,2)
    clf
    ax1 = subplot(1,2,1);
    hold on
    
    plotData = detrend(data(1,1:sizePLI*500));
    timeWindowIdxBoundsSec = timeWindowIdxBounds./60-5;
    patchx = [timeWindowIdxBoundsSec(1,tmpFrame), timeWindowIdxBoundsSec(2,tmpFrame), timeWindowIdxBoundsSec(2,tmpFrame), timeWindowIdxBoundsSec(1,tmpFrame)];
    patchy = [2*min(plotData), 2*min(plotData), 2*max(plotData), 2*max(plotData)];
    
    pData = patch(patchx, patchy,1);
    
    set(pData, 'EdgeColor', 'none')
    set(pData, 'FaceColor', [0.5,0.5,0.5])
    set(pData, 'FaceAlpha', 0.25)
    
    
    
    plot(linspace(-5,5,sizePLI*500), plotData)
    title('Detrended Data from Electrode 1')
    xlabel('Time, minutes')
    ylabel('Voltage, uV')
    ylim([1.1*min(plotData), 1.1*max(plotData)])
    box on
    
    ax2 = subplot(1,2,2);
    
    hold on
    mapCol = gray(128);
    border = 0;
    
    layout = reshape(1:64,8,8);
    
    
    for ii = 1:4
        
        % Find Lower Left corner
        chans = 1:64;
        [~,chanLLIdx] = intersect(layout,chans);
        
        % Subdivide each section into [x,y] sections. Where x and y are layout
        % dimensions.
        rowsFix = linspace(0+border,1-border,size(layout,1)+1);
        
        colsFix = linspace(0+border,1-border,size(layout,2)+1);
        
        for jj = 1:length(chanLLIdx)
            
            chanIdx = layout(chanLLIdx(jj));
            [offx, offy] = ind2sub([size(layout,1),size(layout,2)],find(layout==chanIdx));
            
            xSubPos = [rowsFix(1,1), rowsFix(1,end), rowsFix(1,end), rowsFix(1,1)  ] + offx;
            ySubPos = [colsFix(1,1), colsFix(1,1)  , colsFix(1,end), colsFix(1,end)] + offy;
            
            chanIdx2 = layout(chanLLIdx(jj));
            
            [xpos, ypos] = ind2sub([size(layout,1),size(layout,2)],find(layout==chanIdx2));
            
            colorPatch = mapCol(floor(integralDiff(jj,tmpFrame) * 127+1),:);
            
            pData  = patch( xSubPos, ySubPos, colorPatch);
            
            
            xlim([1,size(layout,1)+1])
            ylim([1,size(layout,2)+1])
            
            set(pData, 'EdgeColor', 'k')
            axis square
            
            set(gca, 'xtick', [1.5:8.5])
            set(gca, 'xticklabel', [1:8])
            
            set(gca, 'ytick', [1.5:8.5])
            set(gca, 'yticklabel', [1:8:57])
            
            xlabel('Channel')
            ylabel('Channel')
            box on
            
        end % END FOR num LL corners
    end    
    
    set(gcf, 'units', 'inches')
    set(gcf, 'Position', [2    2   15    5])
    set(gcf, 'PaperPosition', [2    2   15    5])
    
    pos1=get(ax1, 'Position');
    axes(ax1)
    set(ax1, 'Position', [.055 pos1(2) .575 pos1(4)])
%     
    pos2=get(ax2, 'Position');
    axes(ax2)
    set(ax2, 'Position', [.675 pos2(2) pos2(3) pos2(4)])
    
    drawnow
    
    title('Norm 0.3228')
    
    tmpNumber = num2str(tmpFrame);
    if length(tmpNumber) == 1;
        tmpNumber = ['00', tmpNumber];
    elseif length(tmpNumber) == 2
        tmpNumber = ['0', tmpNumber];
    end
        
    print(['2014PP04Sz4', '_OnsetRef_1to1_5Sec_Frame_', tmpNumber, '.png'],'-dpng')
    
    disp(tmpFrame)
end
