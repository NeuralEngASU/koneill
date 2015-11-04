%% PLI Direction Plot

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
%% 3x3 PLI direction

timeWindowBounds = [[-5:0.25:4.75]; [-4.75:0.25:5]];
time = linspace(-5,5,600);
for tw = 1:size(timeWindowBounds,2)
    
    timeWindow = time >= timeWindowBounds(1,tw) & time < timeWindowBounds(2,tw);
    timeSize = 60* diff(timeWindowBounds(:,tw));
    
    % Channel, time window
end % END FOR each timeWindow
%%

% scaleFac = nan(64,64,600);
% signFac = nan(64,64,600);

u = zeros(64,1,600);
v = zeros(64,1,600);
refChan = 65;

for desChan = 1:64
    
    % Find Channels to Plot (3x3 grid around chanIdx)
    gridDim = [8,8];
    [gridI, gridJ] = ind2sub(gridDim, desChan);
    
    iiNew = [gridI-1, gridI, gridI+1];
    jjNew = [gridJ-1, gridJ, gridJ+1];
    
    % Find index pairs
    [idx2, idx1] = find(true(numel(iiNew),numel(jjNew)));
    indPairs = [reshape(iiNew(idx1), [], 1), reshape(jjNew(idx2), [], 1)];
    
    % Remove invalid pairs (where there is a 0 index)
    indParisIdx = (indPairs(:,1) > 0) & (indPairs(:,2) > 0) & (indPairs(:,1) <= gridDim(1)) & (indPairs(:,2) <= gridDim(2));
    indPairs2 = indPairs(indParisIdx,:);
    
    % Convert index notation to channel number
    chansPlot = sort(sub2ind([8,8], indPairs2(:,1), indPairs2(:,2)), 'ascend');
    % chansPlot = 1:64;
%     desiredRef = [desChan];
    desiredRef = refChan;
    
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
    
    % Grid layout
    layout = reshape(1:64,8,8);
    
    chansPlot = chansPlot(chansPlot~=desChan);
    
%     [x,y] = ind2sub(size(layout),1:64); % Generate x,y coordinates
    [x,y] = ind2sub(size(layout),chansPlot); % Generate x,y coordinates
    [adj(1), adj(2)] = ind2sub(size(layout),desChan);
    
    vectorAngles = atan2(y-adj(2),x-adj(1));% *(180/pi); % Generate vector angles at those coordinates
        
    
    vectorAngles = repmat(vectorAngles',1,1,600);
    
%     twIdx = tw;%(tw-1)*15+1 : (tw-1)*15+15;
    
    scaleFac = [];
    scaleFac(desChan,:,:) = permute(szPLI(:,idx), [3,2,1]);
    %     scaleFacTmp = nan(1,1,64);
    %     scaleFacTmp(:,1,chansPlot(chansPlot~=desChan)) = mean(permute(szPLI(twIdx,idx),[3,1,2]),2);
    
    phiPairs = [repmat(desChan,length(chansPlot),1), chansPlot];
    
    signFac = [];
    
    for ii = 1:length(chansPlot)
        deltaPhi = szPhi(:,:,phiPairs(ii,2)) - szPhi(:,:,phiPairs(ii,1));
        
        % Shift points so they lie within -pi:pi
        deltaPhi(deltaPhi < -pi) = deltaPhi(deltaPhi < -pi) + 2*pi;
        deltaPhi(deltaPhi >  pi) = deltaPhi(deltaPhi >  pi) - 2*pi;
        
        deltaPhi = permute(deltaPhi,[1,3,2]);
        
        signFac(desChan,ii,:) = sign(mean(deltaPhi));
        
    end % END FOR channel
    
    u(desChan,:,:) = mean(cos(vectorAngles).*scaleFac(desChan,:,:) .* signFac(desChan,:,:),2);
    v(desChan,:,:) = mean(sin(vectorAngles).*scaleFac(desChan,:,:) .* signFac(desChan,:,:),2);
    
    disp(desChan)
    
end % END desChan

%% Integrate

winSize = 5;
winSizeAdj = winSize-1;

timeWindowIdxBounds = [[240:600-winSizeAdj];[240+winSizeAdj:600]];

integralDiff = [];

% timeWindowBounds = [[-5:0.25:4.75]; [-4.75:0.25:5]];

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
% desiredRef = [chanIdx];
desiredRef = refChan;

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


%% Break into windows



% for winNum = 1:size(timeWindowIdxBounds,2)
%     
%     twIdx = (winNum-1)*15+1 : (winNum-1)*15+15;
%     
%     uWin(:,:,winNum) = nanmean(u(:,:,twIdx),3);
%     vWin(:,:,winNum) = nanmean(v(:,:,twIdx),3);
% end % END FOR

for winNum = 1:size(timeWindowIdxBounds,2)

    twIdx = timeWindowIdxBounds(1,winNum):timeWindowIdxBounds(2,winNum);

    uWin(:,:,winNum) = nanmean(u(:,:,twIdx),3);
    vWin(:,:,winNum) = nanmean(v(:,:,twIdx),3);
end % END FOR



%% Find Max

uWinUnwrap = uWin(:);
vWinUnwrap = vWin(:);

dist = sqrt(uWinUnwrap.^2 + vWinUnwrap.^2);
[maxDist, maxDistIdx] = max(dist)

uWinScale = uWin./max(dist);
vWinScale = vWin./max(dist);

%% Animate

for tmpFrame = 1:size(uWinScale,3)
    
    close all
    figure(1);
    hold on
    mapCol = gray(128);
    border = 0;
    
    layout = reshape(1:64,8,8);
    
    [x,y] = ind2sub(size(layout),1:64);
    
    x = x';
    y = y';
    
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
            
            colorPatch = mapCol(floor(128),:);
            
            pData  = patch( xSubPos, ySubPos, colorPatch);
            
            
            xlim([0.5,size(layout,1)+1.5])
            ylim([0.5,size(layout,2)+1.5])
            
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
    
    xVals = x - uWinScale(:,1,tmpFrame)/2 + 0.5;
    yVals = y - vWinScale(:,1,tmpFrame)/2 + 0.5;
    
    uVals = uWinScale(:,1,tmpFrame);
    vVals = vWinScale(:,1,tmpFrame);
    
    quiver(xVals, yVals, uVals, vVals,0);
    
    print(['2014PP04Sz5', '_Direction_5Sec_Frame_', num2str(tmpFrame), '.png'],'-dpng')
    
    disp(tmpFrame)
end

%% Animate and Color (Use SfN_PosterPlots line ~625)

for tmpFrame = 1:size(uWinScale,3)
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
    
    [x,y] = ind2sub(size(layout),1:64);
    
    x = x';
    y = y';
    
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
    
    xVals = x - uWinScale(:,1,tmpFrame)/2 + 0.5;
    yVals = y - vWinScale(:,1,tmpFrame)/2 + 0.5;
    
    uVals = uWinScale(:,1,tmpFrame);
    vVals = vWinScale(:,1,tmpFrame);
    
    qh = quiver(xVals, yVals, uVals, vVals,0, 'color', 'w');
    
    set(qh,'linewidth',1);
    
    
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
    
    tmpNumber = num2str(tmpFrame);
    if length(tmpNumber) == 1;
        tmpNumber = ['00', tmpNumber];
    elseif length(tmpNumber) == 2
        tmpNumber = ['0', tmpNumber];
    end
        
    print(['2012PP05Sz1', '_Direction_5Sec_Frame_', tmpNumber, '.png'],'-dpng')
    
    disp(tmpFrame)
end

%% Append to AVI

writerObj = VideoWriter('2014PP04Sz4_OnsetRef_1to1_5SecBins.avi');
writerObj.FrameRate=10;
open(writerObj);



expr = '(2012PP05Sz1_5Sec_Frame_[0-9]+\.png)';

imFileNames = dir(['D:\PLI\SeizureDetection\AnimatedFigures\2014PP04Sz4_OnsetRef\OnsetRef_1to1\*.png']);

for ii = 1:117
        
    
    imData = imread(imFileNames(ii).name);
    writeVideo(writerObj, imData);
    
end % END FOR

close(writerObj);

%% GIF

for ii = 1:40
        
   
    
    imData = imread(imFileNames(ii).name);
    
    [A,map] = rgb2ind(imData,256);
    if ii == 1;
        imwrite(A,map,'2014PP04Sz5_Direction_15SecBins.gif','gif','LoopCount',Inf,'DelayTime',0.5);
    else
        imwrite(A,map,'2014PP04Sz5_Direction_15SecBins.gif','gif','WriteMode','append','DelayTime',0.5);
    end
    
end % END FOR


% EOF