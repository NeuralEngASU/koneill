%%
szFilePath = 'D:\PLI\SeizureDetection\Sz\CAR';
nonSzFilePath = 'D:\PLI\SeizureDetection\NonSz\CAR';

%% Plot average PLI over 64 chans

szFilePath = 'D:\PLI\SeizureDetection\Sz\CAR';
nonSzFilePath = 'D:\PLI\SeizureDetection\NonSz\CAR';

szFileName = dir([szFilePath, '\*.mat']);
nonSzFileName = dir([nonSzFilePath, '\*.mat']);
numFiles = size(szFileName,1);

fileOfInterest = 1;

for fileOfInterest=1:1%68
%%
% Load Sz
load(fullfile(szFilePath, szFileName(fileOfInterest).name));
szChanPairNums = chanPairNums;
szHeader = Header;
szPLI = p;
szR = r;
szParams = params;
szPhi = phi;

% Load NonSz
load(fullfile(nonSzFilePath, nonSzFileName(fileOfInterest).name));
nonSzChanPairNums = chanPairNums;
nonSzHeader = Header;
nonSzPLI = p;
nonSzR = r;
nonSzParams = params;
nonSzPhi = phi;

% Clear old data
clear chanPairNums Header p params phi r

%

chansPlot = [1:64];
desiredRef = [];

% Define channel pairs

if isempty(desiredRef)
     desiredChanPairs = nchoosek(sort(unique(chansPlot),'ascend'),2);
else
    desiredChanPairs = [repmat(desiredRef, size(chansPlot(:),1), 1), chansPlot(:)];
end % END IF isempty(desiredRef)

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

% Calculate Error
gridErrSz = std(szPLI(1:sizePLI,idx), 0, 2);
gridErrSz = gridErrSz/sqrt(length(idx)); % Standard error
gridErrSz = 2*gridErrSz; % 2*standard error.

gridErrNonSz = std(nonSzPLI(1:sizePLI,idx), 0, 2);
gridErrNonSz = gridErrNonSz/sqrt(length(idx)); % Standard error
gridErrNonSz = 2*gridErrNonSz; % 2*standard error.

timeMax = size(szPLI,1);
sizeMax = szHeader.params.Fs*timeMax;

patchx = linspace(-5,5,sizePLI);
patchxx = [patchx, fliplr(patchx)];

% Plot
figure;
subplot(2,1,1) % Seizure Plot
patchdata =  [[meanSzPLI' + gridErrSz]', fliplr([meanSzPLI' - gridErrSz]')];

pData = patch(patchxx, patchdata, 1);
hold on
lData = plot(patchx, meanSzPLI, 'b', 'linewidth', 1);
hold off
set(pData, 'FaceColor', 'k')
set(pData, 'EdgeColor', 'none')
set(pData, 'FaceAlpha', 0.25)
set(gca, 'XTick', unique([0:10:patchx(end), patchx(end)]))

ylim([0,1])
t = title(['Sz: ', szFileName(fileOfInterest).name]);
ylabel('WPLI')
xlabel('Time, min')

set(t,'Interpreter','none');

% Subplot 2
subplot(2,1,2) % Non Seizure plot
patchdata =  [[meanNonSzPLI' + gridErrNonSz]', fliplr([meanNonSzPLI' - gridErrNonSz]')];

pData = patch(patchxx, patchdata, 1);
hold on
lData = plot(patchx, meanNonSzPLI, 'b', 'linewidth', 1);
hold off
set(pData, 'FaceColor', 'k')
set(pData, 'EdgeColor', 'none')
set(pData, 'FaceAlpha', 0.25)
set(gca, 'XTick', unique([0:10:patchx(end), patchx(end)]))

ylim([0,1])
t = title(['NonSz: ', nonSzFileName(fileOfInterest).name]);
ylabel('WPLI')
xlabel('Time, min')

set(t,'Interpreter','none');

print([szFileName(fileOfInterest).name, '_WholeGrid.png'],'-dpng')
% close all

end

%% FOR loop over grid. Compare 3x3 square.

szFilePath = 'D:\PLI\SeizureDetection\Sz\HilbertFirst';
nonSzFilePath = 'D:\PLI\SeizureDetection\NonSz\HilbertFirst';

szFileName = dir([szFilePath, '\*.mat']);
nonSzFileName = dir([nonSzFilePath, '\*.mat']);
numFiles = size(szFileName,1);

fileOfInterest = 58;

% Load Sz
% load(fullfile(szFilePath, szFileName(fileOfInterest).name));
szFileName(fileOfInterest).name = '2012PP05Sz1_PLI_winSize1.mat';
load(fullfile(szFilePath, '2012PP05Sz1_PLI_winSize1.mat'));
szChanPairNums = chanPairNums;
szHeader = Header;
szPLI = p;
szR = r;
szParams = params;
szPhi = phi;

% Load NonSz
% load(fullfile(nonSzFilePath, nonSzFileName(fileOfInterest).name));
nonSzFileName(fileOfInterest).name = '2012PP05NonSz1_PLI_winSize1.mat';
load(fullfile(nonSzFilePath, '2012PP05NonSz1_PLI_winSize1.mat'));

nonSzChanPairNums = chanPairNums;
nonSzHeader = Header;
nonSzPLI = p;
nonSzR = r;
nonSzParams = params;
nonSzPhi = phi;

% Clear old data
clear chanPairNums Header p params phi r

for chanIdx = 10

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

% Calculate Error
gridErrSz = std(szPLI(1:sizePLI,idx), 0, 2);
gridErrSz = gridErrSz/sqrt(length(idx)); % Standard error
gridErrSz = 2*gridErrSz; % 2*standard error.

gridErrNonSz = std(nonSzPLI(1:sizePLI,idx), 0, 2);
gridErrNonSz = gridErrNonSz/sqrt(length(idx)); % Standard error
gridErrNonSz = 2*gridErrNonSz; % 2*standard error.

timeMax = size(szPLI,1);
sizeMax = szHeader.params.Fs*timeMax;

patchx = linspace(-5,5,sizePLI);
patchxx = [patchx, fliplr(patchx)];

% Plot
figure;
subplot(2,1,1) % Seizure Plot
patchdata =  [[meanSzPLI' + gridErrSz]', fliplr([meanSzPLI' - gridErrSz]')];

pData = patch(patchxx, patchdata, 1);
hold on
lData = plot(patchx, meanSzPLI, 'b', 'linewidth', 1);
hold off
set(pData, 'FaceColor', 'k')
set(pData, 'EdgeColor', 'none')
set(pData, 'FaceAlpha', 0.25)
set(gca, 'XTick', unique([0:10:patchx(end), patchx(end)]))

ylim([0,1])
% t = title(['Sz: ', szFileName(fileOfInterest).name, '   ', num2str(chansPlot')]);
title('2012PP05 Seizure Clip 1. Channel 10', 'FontSize', 25)
ylabel('WPLI', 'FontSize', 25)
xlabel('Time, min', 'FontSize', 25)

set(gca, 'xtick', [-5:5])
set(gca, 'xticklabels', [-5:5], 'FontSize', 20)
box on

% set(t,'Interpreter','none');

% Subplot 2
subplot(2,1,2) % Non Seizure plot
patchdata =  [[meanNonSzPLI' + gridErrNonSz]', fliplr([meanNonSzPLI' - gridErrNonSz]')];

pData = patch(patchxx, patchdata, 1);
hold on
lData = plot(patchx, meanNonSzPLI, 'b', 'linewidth', 1);
hold off
set(pData, 'FaceColor', 'k')
set(pData, 'EdgeColor', 'none')
set(pData, 'FaceAlpha', 0.25)
set(gca, 'XTick', unique([0:10:patchx(end), patchx(end)]))

ylim([0,1])
xlim([-5,5])
% t = title(['NonSz: ', nonSzFileName(fileOfInterest).name]);
title('2012PP05 Non-Seizure Clip 1. Channel 10', 'FontSize', 25)
ylabel('WPLI', 'FontSize', 25)
xlabel('Time, min', 'FontSize', 25)
set(gca, 'xtick', [-5:5])
set(gca, 'xticklabels', [-5:5], 'FontSize', 20)
% set(t,'Interpreter','none');

% print([szFileName(fileOfInterest).name, '_3x3Grid_Chan', num2str(chanIdx), '.png'],'-dpng')
% close all

    set(gcf, 'units', 'inches')
    set(gcf, 'Position', [2    2   16    8])
    set(gcf, 'PaperPosition', [2    2   16    8])
box on
end

%% FOR Loop, compare grid to everything else

%% FOR loop over grid. Compare full grid.

szFilePath = 'D:\PLI\SeizureDetection\Sz\HilbertFirst';
nonSzFilePath = 'D:\PLI\SeizureDetection\NonSz\HilbertFirst';

szFileName = dir([szFilePath, '\*.mat']);
nonSzFileName = dir([nonSzFilePath, '\*.mat']);
numFiles = size(szFileName,1);

fileOfInterest = 1;

% Load Sz
% load(fullfile(szFilePath, szFileName(fileOfInterest).name));
szFileName(fileOfInterest).name = '2014PP06Sz6_PLI_winSize1.mat';
load(fullfile(szFilePath, '2014PP06Sz6_PLI_winSize1.mat'));
szChanPairNums = chanPairNums;
szHeader = Header;
szPLI = p;
szR = r;
szParams = params;
szPhi = phi;

% Load NonSz
% load(fullfile(nonSzFilePath, nonSzFileName(fileOfInterest).name));
nonSzFileName(fileOfInterest).name = '2014PP06NonSz6_PLI_winSize1.mat';
load(fullfile(nonSzFilePath, '2014PP06NonSz6_PLI_winSize1.mat'));

nonSzChanPairNums = chanPairNums;
nonSzHeader = Header;
nonSzPLI = p;
nonSzR = r;
nonSzParams = params;
nonSzPhi = phi;

% Clear old data
clear chanPairNums Header p params phi r

for chanIdx = 1:64

chansPlot = [1:64]
desiredRef = [];

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



for ii = 1:length(idx)
    tmpIdx = idx(ii);
    smoothSzPLI(:,ii) = smooth(szPLI(:,tmpIdx));
    smoothNonSzPLI(:,ii) = smooth(nonSzPLI(:,tmpIdx));
end % END FOR

meanSzPLI = mean(smoothSzPLI,2);
meanNonSzPLI = mean(smoothNonSzPLI,2);

% Calculate Error
gridErrSz = std(smoothSzPLI, 0, 2);
gridErrSz = gridErrSz/sqrt(length(idx)); % Standard error
gridErrSz = 2*gridErrSz; % 2*standard error.

gridErrNonSz = std(smoothNonSzPLI, 0, 2);
gridErrNonSz = gridErrNonSz/sqrt(length(idx)); % Standard error
gridErrNonSz = 2*gridErrNonSz; % 2*standard error.

timeMax = size(szPLI,1);
sizeMax = szHeader.params.Fs*timeMax;

patchx = linspace(-5,5,sizePLI);
patchxx = [patchx, fliplr(patchx)];

% Plot
figure;
subplot(2,1,1) % Seizure Plot
patchdata =  [[meanSzPLI + gridErrSz]', fliplr([meanSzPLI - gridErrSz]')];

pData = patch(patchxx, patchdata, 1);
hold on
lData = plot(patchx, meanSzPLI, 'b', 'linewidth', 1);
hold off
set(pData, 'FaceColor', 'k')
set(pData, 'EdgeColor', 'none')
set(pData, 'FaceAlpha', 0.25)
set(gca, 'XTick', unique([0:10:patchx(end), patchx(end)]))

ylim([0,1])
t = title(['Sz: ', szFileName(fileOfInterest).name, ' Chan: ', num2str(desiredRef)]);
ylabel('WPLI')
xlabel('Time, min')

set(t,'Interpreter','none');

% Subplot 2
subplot(2,1,2) % Non Seizure plot
patchdata =  [[meanNonSzPLI + gridErrNonSz]', fliplr([meanNonSzPLI - gridErrNonSz]')];

pData = patch(patchxx, patchdata, 1);
hold on
lData = plot(patchx, meanNonSzPLI, 'b', 'linewidth', 1);
hold off
set(pData, 'FaceColor', 'k')
set(pData, 'EdgeColor', 'none')
set(pData, 'FaceAlpha', 0.25)
set(gca, 'XTick', unique([0:10:patchx(end), patchx(end)]))

ylim([0,1])
t = title(['NonSz: ', nonSzFileName(fileOfInterest).name]);
ylabel('WPLI')
xlabel('Time, min')

set(t,'Interpreter','none');

print([szFileName(fileOfInterest).name, '_FullGrid_Smooth_Chan', num2str(chanIdx), '.png'],'-dpng')
close all

end

%% Integral Plots of 3x3 subgrids

szFilePath = 'D:\PLI\SeizureDetection\Sz\HilbertFirst';
nonSzFilePath = 'D:\PLI\SeizureDetection\NonSz\HilbertFirst';

szFileName = dir([szFilePath, '\*.mat']);
nonSzFileName = dir([nonSzFilePath, '\*.mat']);
numFiles = size(szFileName,1);

fileOfInterest = 1;


szFileName(fileOfInterest).name = '2012PP05Sz1_PLI_winSize1.mat';
nonSzFileName(fileOfInterest).name = '2012PP05NonSz1_PLI_winSize1.mat';

% Load Sz
load(fullfile(szFilePath, szFileName(fileOfInterest).name));
szChanPairNums = chanPairNums;
szHeader = Header;
szPLI = p;
szR = r;
szParams = params;
szPhi = phi;

% Load NonSz
load(fullfile(nonSzFilePath, nonSzFileName(fileOfInterest).name));
nonSzChanPairNums = chanPairNums;
nonSzHeader = Header;
nonSzPLI = p;
nonSzR = r;
nonSzParams = params;
nonSzPhi = phi;

% Clear old data
clear chanPairNums Header p params phi r
%%
integralDiff = [];
badChan = [33,32]; % For 2012PP05Sz1
chanLoop = 1:64;

for ii = 1:length(badChan)
    chanLoop = chanLoop(chanLoop ~= badChan(ii));
end % END FOR

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

meanSzPLI = mean(szPLI(1:sizePLI,idx)');
meanNonSzPLI = mean(nonSzPLI(1:sizePLI,idx)');

timeMax = size(szPLI,1);
sizeMax = szHeader.params.Fs*timeMax;

time = linspace(-5,5,sizePLI);

% Channel, time window, sz/nonSz
rawIntegral(chanIdx,1,1) = sum(meanSzPLI)/sizePLI;
rawIntegral(chanIdx,1,2) = sum(meanNonSzPLI)/sizePLI;

% timeWindowBounds = [[-2.33;-2], [0.5;0.88], [1.33;1.66], [3.66; 4]]; % 2014PP04Sz4
timeWindowBounds = [[-2.33;-2], [0.33;0.66], [1.66;2], [3.66; 4]]; % 2012PP05Sz1


timeWindow(:,1) = time >= timeWindowBounds(1,1) & time < timeWindowBounds(2,1);
timeWindow(:,2) = time >= timeWindowBounds(1,2) & time < timeWindowBounds(2,2);
timeWindow(:,3) = time >= timeWindowBounds(1,3) & time < timeWindowBounds(2,3);
timeWindow(:,4) = time >= timeWindowBounds(1,4) & time < timeWindowBounds(2,4);



% Channel, time window
% integralDiff(chanIdx,1) = abs(sum(meanSzPLI(timeWindow(:,1)))/60 - sum(meanNonSzPLI(timeWindow(:,1)))/60);
% integralDiff(chanIdx,2) = abs(sum(meanSzPLI(timeWindow(:,2)))/60 - sum(meanNonSzPLI(timeWindow(:,2)))/60);
% integralDiff(chanIdx,3) = abs(sum(meanSzPLI(timeWindow(:,3)))/60 - sum(meanNonSzPLI(timeWindow(:,3)))/60);
% integralDiff(chanIdx,4) = abs(sum(meanSzPLI(timeWindow(:,4)))/60 - sum(meanNonSzPLI(timeWindow(:,4)))/60);

% Channel, time window
integralDiff(chanIdx,1) = abs(sum(meanSzPLI(timeWindow(:,1)))/60);
integralDiff(chanIdx,2) = abs(sum(meanSzPLI(timeWindow(:,2)))/60);
integralDiff(chanIdx,3) = abs(sum(meanSzPLI(timeWindow(:,3)))/60);
integralDiff(chanIdx,4) = abs(sum(meanSzPLI(timeWindow(:,4)))/60);


end % END FOR

% integralDiff = permute(integralDiff, [1,3,2]);

% integralDiff = reshape(integralDiff, 8,8,4);
% integralDiff = integralDiff ./ max(integralDiff(:));
%


% maxVal = 0.3211; % 2012PP05Sz1
maxVal = 0.33;
% maxVal = 0.1937; % 2014PP04Sz4
integralDiff = integralDiff ./ maxVal;

load('E:\data\human CNS\EMD\Sz\clips\2012PP05Sz1.mat')

% integralDiff2 = integralDiff;
% integralDiff = integralDiff ./ max(integralDiff(:));

figure;
mapCol = gray(128);
border = 0;

layout = reshape(1:64,8,8);

plotData = detrend(data(10,1:sizePLI*500))./1000;

subplot(2,1,1)
hold on
box on
for kk = 1:4
    
    patchx = [timeWindowBounds(1,kk), timeWindowBounds(2,kk), timeWindowBounds(2,kk), timeWindowBounds(1,kk)];
    patchy = [2*min(plotData), 2*min(plotData), 2*max(plotData), 2*max(plotData)];
    
    pData = patch(patchx, patchy,1);
    
    set(pData, 'EdgeColor', 'none')
    set(pData, 'FaceColor', [0.5,0.5,0.5])
    set(pData, 'FaceAlpha', 0.25)
    
end % END FOR


plot(linspace(-5,5,sizePLI*500), plotData)
title('Detrended Data from Electrode 10', 'FontSize', 25)
xlabel('Time, minutes', 'FontSize', 25)
ylabel('Voltage, mV', 'FontSize', 25)
ylim([1.1*min(plotData), 1.1*max(plotData)])

a = get(gca, 'xticklabels');
b = get(gca, 'yticklabels');
set(gca, 'xticklabels', a, 'FontSize', 20)
set(gca, 'yticklabels', b, 'FontSize', 20)



for ii = 1:4
   
    subplot(2,4,ii+4)
    
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
        
        colorPatch = mapCol(floor(integralDiff(jj,ii) * 127+1),:);
        
        pData  = patch( xSubPos, ySubPos, colorPatch);
        
        
        xlim([1,size(layout,1)+1])
        ylim([1,size(layout,2)+1])
        
        set(pData, 'EdgeColor', 'k')
        axis square
        
        titleString = {'Pre-Ictal', 'Ictal', 'Early Post-Ictal', 'Late Post-Ictal'};
        tmpTitle = sprintf([titleString{ii}, '\n', num2str(timeWindowBounds(1,ii)), ' to ', num2str(timeWindowBounds(2,ii)), ' minutes']);
        title(tmpTitle, 'FontSize', 25) 
        
        set(gca, 'xtick', [1.5:8.5])
        set(gca, 'xticklabel', [1:8], 'FontSize', 20)
        
        set(gca, 'ytick', [1.5:8.5])
        set(gca, 'yticklabel', [1:8:57], 'FontSize', 20)
        
        xlabel('Channel')
        ylabel('Channel')
        box on
        
    end % END FOR num LL corners    
end

set(gcf, 'units', 'inches')
set(gcf, 'Position', [2    2   16    8])
set(gcf, 'PaperPosition', [2    2   16    8])

%% Animated Integral Plots of 3x3 Grid 

szFilePath = 'D:\PLI\SeizureDetection\Sz\HilbertFirst';
nonSzFilePath = 'D:\PLI\SeizureDetection\NonSz\HilbertFirst';

szFileName = dir([szFilePath, '\*.mat']);
nonSzFileName = dir([nonSzFilePath, '\*.mat']);
numFiles = size(szFileName,1);

fileOfInterest = 1;

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
    integralDiff(chanIdx,tw) = abs(sum(meanSzPLI(timeWindow))/timeSize - sum(meanNonSzPLI(timeWindow))/timeSize);
    
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

%% Append Gif

expr = '(2014PP04Sz4_15Sec_Frame_[0-9]+\.png)';

imFileNames = dir(['D:\PLI\SeizureDetection\AnimatedFigures', '\*.png']);

for ii = 1:40
    listNum = 61:100;
    
    imIdx = listNum(ii);
    
    imData = imread(imFileNames(imIdx).name);
    
    [A,map] = rgb2ind(imData,256);
    if ii == 1;
        imwrite(A,map,'2014PP04Sz4_15Sec_gif.gif','gif','LoopCount',Inf,'DelayTime',0.5);
    else
        imwrite(A,map,'2014PP04Sz4_15Sec_gif.gif','gif','WriteMode','append','DelayTime',0.5);
    end
    
end % END FOR

%% PLI Direction Plot

szFilePath = 'D:\PLI\SeizureDetection\Sz\HilbertFirst';
nonSzFilePath = 'D:\PLI\SeizureDetection\NonSz\HilbertFirst';

szFileName = dir([szFilePath, '\*.mat']);
nonSzFileName = dir([nonSzFilePath, '\*.mat']);
numFiles = size(szFileName,1);

fileOfInterest = 1;

% Load Sz
% load(fullfile(szFilePath, szFileName(fileOfInterest).name));
szFileName(fileOfInterest).name = '2014PP04Sz7_PLI_winSize1.mat';
load(fullfile(szFilePath, szFileName(fileOfInterest).name));
szChanPairNums = chanPairNums;
szHeader = Header;
szPLI = p;
szR = r;
szParams = params;
szPhi = phi;

% Load NonSz
% load(fullfile(nonSzFilePath, nonSzFileName(fileOfInterest).name));
nonSzFileName(fileOfInterest).name = '2014PP04NonSz7_PLI_winSize1.mat';
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
load('E:\data\human CNS\EMD\Sz\clips\2014PP04Sz4.mat')

%%

% Calculate deltaPhi
phi1 = szPhi(:,:,1);
phi2 = szPhi(:,:,2);

deltaPhi = phi1 - phi2;

% Shift points so they lie within -pi:pi
deltaPhi(deltaPhi < -pi) = deltaPhi(deltaPhi < -pi) + 2*pi;
deltaPhi(deltaPhi >  pi) = deltaPhi(deltaPhi >  pi) - 2*pi;

% deltaPhi = deltaPhi(:);

t = linspace(-5,5,length(deltaPhi));

plot(t,deltaPhi)

%%
% Subset PLI

for desChan = 1:64
    
    %     % Find Channels to Plot (3x3 grid around chanIdx)
    %     gridDim = [8,8];
    %     [gridI, gridJ] = ind2sub(gridDim, desChan);
    %
    %     iiNew = [gridI-1, gridI, gridI+1];
    %     jjNew = [gridJ-1, gridJ, gridJ+1];
    %
    %     % Find index pairs
    %     [idx2, idx1] = find(true(numel(iiNew),numel(jjNew)));
    %     indPairs = [reshape(iiNew(idx1), [], 1), reshape(jjNew(idx2), [], 1)];
    %
    %     % Remove invalid pairs (where there is a 0 index)
    %     indParisIdx = (indPairs(:,1) > 0) & (indPairs(:,2) > 0) & (indPairs(:,1) <= gridDim(1)) & (indPairs(:,2) <= gridDim(2));
    %     indPairs2 = indPairs(indParisIdx,:);
    %
    %     % Convert index notation to channel number
    %     chansPlot = sort(sub2ind([8,8], indPairs2(:,1), indPairs2(:,2)), 'ascend');
    chansPlot = 1:64;
    desiredRef = [desChan];
    
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
    
    layout = reshape(1:64,8,8);
    
    [x,y] = ind2sub(size(layout),1:64); % Generate x,y coordinates
    [adj(1), adj(2)] = ind2sub(size(layout),desChan);
    
    vectorAngles = atan2(y-adj(2),x-adj(1));% *(180/pi); % Generate vector angles at those coordinates
    
    [~,maxIdx] = max(szPLI(:,1));
    
    maxVal = 1;
    
    scaleFacTmp = szPLI(maxIdx, idx);
    
    %     scaleFacTmp = nan(1,64);
    %     scaleFacTmp(chansPlot) = szPLI(maxIdx,idx);
    
    
    if desChan == 1
        scaleFacTmp = [NaN, scaleFacTmp];
    elseif desChan == 64
        scaleFacTmp = [scaleFacTmp, NaN];
    else
        scaleFacTmp = [scaleFacTmp(1:(desChan-1)), NaN, scaleFacTmp(desChan:end)];
    end % END IF
    
    
    scaleFac(desChan,:) = scaleFacTmp;
    
    phiPairs = [repmat(desChan,64,1), [1:64]'];
    
    for ii = 1:64
        deltaPhi = szPhi(:,maxIdx,phiPairs(ii,1)) - szPhi(:,maxIdx,phiPairs(ii,2));
        
        % Shift points so they lie within -pi:pi
        deltaPhi(deltaPhi < -pi) = deltaPhi(deltaPhi < -pi) + 2*pi;
        deltaPhi(deltaPhi >  pi) = deltaPhi(deltaPhi >  pi) - 2*pi;
        
        signFac(desChan,ii) = sign(mean(deltaPhi));
        
    end
    
    u(desChan,:) = cos(vectorAngles).*scaleFac(desChan,:) .* signFac(desChan,:);
    v(desChan,:) = sin(vectorAngles).*scaleFac(desChan,:) .* signFac(desChan,:);
    
    
    
end

figure;
quiver(x,y,nanmean(u),nanmean(v));
title(['Coupling Direction Map'])
axis square

figure;
quiver(x,y,u(57,:),v(57,:));
title(['Coupling Direction Map'])
axis square
    
%% 3x3 PLI direction


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
    desiredRef = [desChan];
    
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
    
    layout = reshape(1:64,8,8);
    
    [x,y] = ind2sub(size(layout),1:64); % Generate x,y coordinates
    [adj(1), adj(2)] = ind2sub(size(layout),desChan);

    vectorAngles = atan2(y-adj(2),x-adj(1));% *(180/pi); % Generate vector angles at those coordinates
    
    [~,maxIdx] = max(szPLI(:,1));
    
    maxIdx = 347;
    
    maxVal = 1;
        
%     scaleFacTmp = szPLI(maxIdx, idx);
    
    scaleFacTmp = nan(1,64);
    scaleFacTmp(chansPlot(chansPlot~=desChan)) = szPLI(maxIdx,idx);
    
    
%     if desChan == 1
%         scaleFacTmp = [NaN, scaleFacTmp];
%     elseif desChan == 64
%         scaleFacTmp = [scaleFacTmp, NaN];
%     else
%         scaleFacTmp = [scaleFacTmp(1:(desChan-1)), NaN, scaleFacTmp(desChan:end)];
%     end % END IF
    
    
    scaleFac(desChan,:) = scaleFacTmp;
    
    phiPairs = [repmat(desChan,64,1), [1:64]'];
    
    for ii = 1:64
        deltaPhi = szPhi(:,maxIdx,phiPairs(ii,2)) - szPhi(:,maxIdx,phiPairs(ii,1));
        
        % Shift points so they lie within -pi:pi
        deltaPhi(deltaPhi < -pi) = deltaPhi(deltaPhi < -pi) + 2*pi;
        deltaPhi(deltaPhi >  pi) = deltaPhi(deltaPhi >  pi) - 2*pi;
        
        signFac(desChan,ii) = sign(mean(deltaPhi));
        
    end
    
    u(desChan,:) = cos(vectorAngles).*scaleFac(desChan,:) .* signFac(desChan,:);
    v(desChan,:) = sin(vectorAngles).*scaleFac(desChan,:) .* signFac(desChan,:);
    

    
end

figure;
quiver(x,y,nanmean(u),nanmean(v));
title(['Coupling Direction Map ', num2str(maxIdx)])
axis square
xlim([0,9])
ylim([0,9])

% figure;
% quiver(x,y,u(60,:),v(60,:));
% title(['Coupling Direction Map'])
% axis square

%% 2012PP05, 2014PP04

% fileList{1} = '2012PP05Sz1_PLI_winSize1.mat';
% fileList{2} = '2012PP05Sz3_PLI_winSize1.mat';
% fileList{3} = '2012PP05Sz4_PLI_winSize1.mat';
% fileList{4} = '2012PP05Sz7_PLI_winSize1.mat';
% fileList{5} = '2012PP05Sz8_PLI_winSize1.mat';
% fileList{6} = '2014PP04Sz4_PLI_winSize1.mat';
% fileList{7} = '2014PP04Sz5_PLI_winSize1.mat';
% fileList{8} = '2014PP04Sz6_PLI_winSize1.mat';
% fileList{9} = '2014PP04Sz7_PLI_winSize1.mat';
% fileList{10} = '2014PP04Sz8_PLI_winSize1.mat';
% fileList{11} = '2014PP04Sz9_PLI_winSize1.mat';
% fileList{12} = '2014PP04Sz10_PLI_winSize1.mat';
% fileList{13} = '2014PP04Sz11_PLI_winSize1.mat';
% fileList{14} = '2014PP04Sz12_PLI_winSize1.mat';
% fileList{15} = '2014PP04Sz13_PLI_winSize1.mat';
% fileList{16} = '2014PP04Sz14_PLI_winSize1.mat';
% fileList{17} = '2014PP04Sz15_PLI_winSize1.mat';
% fileList{18} = '2014PP04Sz16_PLI_winSize1.mat';

szFilePath = 'D:\PLI\SeizureDetection\Sz\HilbertFirst';

szFileName = dir([szFilePath, '\*.mat']);
numFiles = size(szFileName,1);
interestingSz = [1,2,3,20,28,30,35,36,37,38,39,40,41,43,69];

% Load Sz
% load(fullfile(szFilePath, szFileName(fileOfInterest).name));
szFileName(fileOfInterest).name = '2014PP04Sz7_PLI_winSize1.mat';


for ii = interestingSz
    
    load(szFileName(ii).name);
    figure(ii)
    plot(mean(p(:,1:2),2))
    title(szFileName(ii).name);
    ylim([0,1])
    
end % END FOR


%% PLI Voltage Plots

szFilePath = 'D:\PLI\SeizureDetection\Sz\HilbertFirst';
nonSzFilePath = 'D:\PLI\SeizureDetection\NonSz\HilbertFirst';

szFileName = dir([szFilePath, '\*.mat']);
nonSzFileName = dir([nonSzFilePath, '\*.mat']);
numFiles = size(szFileName,1);

fileOfInterest = 1;

% Load Sz
% load(fullfile(szFilePath, szFileName(fileOfInterest).name));
szFileName(fileOfInterest).name = '2012PP05Sz1_PLI_winSize1.mat';
load(fullfile(szFilePath, szFileName(fileOfInterest).name));
szChanPairNums = chanPairNums;
szHeader = Header;
szPLI = p;
szR = r;
szParams = params;
szPhi = phi;

% Load NonSz
% load(fullfile(nonSzFilePath, nonSzFileName(fileOfInterest).name));
nonSzFileName(fileOfInterest).name = '2012PP05NonSz1_PLI_winSize1.mat';
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
load('E:\data\human CNS\EMD\Sz\clips\2012PP05Sz1.mat')
szData = data;
load('E:\data\human CNS\EMD\NonSz\clips\2012PP05NonSz1.mat')
nonSzData = data;

%%

time = linspace(-5,5,300000);

% Plot
figure;
subplot(2,1,1) % Seizure Plot
plotData1 = detrend(szData(10,1:300000))./1000;
plot(time, plotData1);

ylim([min(plotData1)*1.1,max(plotData1)*1.1])
title('2012PP05 Seizure Clip 1. Channel 10', 'FontSize', 25)
ylabel('Voltage, mV', 'FontSize', 25)
xlabel('Time, min', 'FontSize', 25)

set(gca, 'xtick', [-5:5])
set(gca, 'xticklabels', [-5:5], 'FontSize',20)
box on

% set(t,'Interpreter','none');

% Subplot 2
plotData2 = detrend(nonSzData(10,1:300000))./1000;
subplot(2,1,2) % Non Seizure plot
plot(time, plotData2);
ylim([min(plotData1)*1.1,max(plotData1)*1.1])
xlim([-5,5])

title('2012PP05 Non-Seizure Clip 1. Channel 10', 'FontSize', 25)
ylabel('Voltage, mV', 'FontSize', 25)
xlabel('Time, min', 'FontSize', 25)
set(gca, 'xtick', [-5:5])
set(gca, 'xticklabels', [-5:5], 'FontSize',20)


    set(gcf, 'units', 'inches')
    set(gcf, 'Position', [2    2   16    8])
    set(gcf, 'PaperPosition', [2    2   16    8])
box on

%% PLI Voltage Plots

szFilePath = 'D:\PLI\SeizureDetection\Sz\HilbertFirst';
nonSzFilePath = 'D:\PLI\SeizureDetection\NonSz\HilbertFirst';

szFileName = dir([szFilePath, '\*.mat']);
nonSzFileName = dir([nonSzFilePath, '\*.mat']);
numFiles = size(szFileName,1);

fileOfInterest = 1;

% Load Sz
% load(fullfile(szFilePath, szFileName(fileOfInterest).name));
szFileName(fileOfInterest).name = '2012PP05Sz1_PLI_winSize1.mat';
load(fullfile(szFilePath, szFileName(fileOfInterest).name));
szChanPairNums = chanPairNums;
szHeader = Header;
szPLI = p;
szR = r;
szParams = params;
szPhi = phi;

% Load NonSz
% load(fullfile(nonSzFilePath, nonSzFileName(fileOfInterest).name));
nonSzFileName(fileOfInterest).name = '2012PP05NonSz1_PLI_winSize1.mat';
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
load('E:\data\human CNS\EMD\Sz\clips\2012PP05Sz1.mat')
szData = data;

%% Color legend

c = linspace(0,1,100);
x = 1:100;
y = 1:100;

colormap gray
scatter(x,y,[],c)

colorbar

cbar_handle = findobj(gcf,'tag','Colorbar')
set(cbar_handle, 'YAxisLocation','Left')
% set(cbar_handle, 'YLabel', 'Average WPLI')
% set(cbar_handle,'fontsize',20);
ylabel(cbar_handle, 'Average WPLI', 'fontsize', 25)
set(cbar_handle, 'YTick',linspace(0,1,5), 'fontsize', 20)



%% 3x3 plot from onset

%% Integral Plots of 3x3 subgrids

szFilePath = 'D:\PLI\SeizureDetection\Sz\HilbertFirst';
nonSzFilePath = 'D:\PLI\SeizureDetection\NonSz\HilbertFirst';

szFileName = dir([szFilePath, '\*.mat']);
nonSzFileName = dir([nonSzFilePath, '\*.mat']);
numFiles = size(szFileName,1);

fileOfInterest = 1;


szFileName(fileOfInterest).name = '2012PP05Sz1_PLI_winSize1.mat';
nonSzFileName(fileOfInterest).name = '2012PP05NonSz1_PLI_winSize1.mat';

% Load Sz
load(fullfile(szFilePath, szFileName(fileOfInterest).name));
szChanPairNums = chanPairNums;
szHeader = Header;
szPLI = p;
szR = r;
szParams = params;
szPhi = phi;

% Load NonSz
load(fullfile(nonSzFilePath, nonSzFileName(fileOfInterest).name));
nonSzChanPairNums = chanPairNums;
nonSzHeader = Header;
nonSzPLI = p;
nonSzR = r;
nonSzParams = params;
nonSzPhi = phi;

% Clear old data
clear chanPairNums Header p params phi r
%%
integralDiff = [];
% badChan = [33,23]; % For 2012PP05Sz1
badChan = []; % For 2014PP04Sz4
chanLoop = 1:64;

for ii = 1:length(badChan)
    chanLoop = chanLoop(chanLoop ~= badChan(ii));
end % END FOR

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

meanSzPLI = mean(szPLI(1:sizePLI,idx)');
meanNonSzPLI = mean(nonSzPLI(1:sizePLI,idx)');

timeMax = size(szPLI,1);
sizeMax = szHeader.params.Fs*timeMax;

time = linspace(-5,5,sizePLI);

% Channel, time window, sz/nonSz
rawIntegral(chanIdx,1,1) = sum(meanSzPLI)/sizePLI;
rawIntegral(chanIdx,1,2) = sum(meanNonSzPLI)/sizePLI;

% timeWindowBounds = [[-2.33;-2], [0.5;0.88], [1.33;1.66], [3.66; 4]]; % 2014PP04Sz4
timeWindowBounds = [[-2.33;-2], [0.33;0.66], [1.66;2], [3.66; 4]]; % 2012PP05Sz1


timeWindow(:,1) = time >= timeWindowBounds(1,1) & time < timeWindowBounds(2,1);
timeWindow(:,2) = time >= timeWindowBounds(1,2) & time < timeWindowBounds(2,2);
timeWindow(:,3) = time >= timeWindowBounds(1,3) & time < timeWindowBounds(2,3);
timeWindow(:,4) = time >= timeWindowBounds(1,4) & time < timeWindowBounds(2,4);



% Channel, time window
% integralDiff(chanIdx,1) = abs(sum(meanSzPLI(timeWindow(:,1)))/60 - sum(meanNonSzPLI(timeWindow(:,1)))/60);
% integralDiff(chanIdx,2) = abs(sum(meanSzPLI(timeWindow(:,2)))/60 - sum(meanNonSzPLI(timeWindow(:,2)))/60);
% integralDiff(chanIdx,3) = abs(sum(meanSzPLI(timeWindow(:,3)))/60 - sum(meanNonSzPLI(timeWindow(:,3)))/60);
% integralDiff(chanIdx,4) = abs(sum(meanSzPLI(timeWindow(:,4)))/60 - sum(meanNonSzPLI(timeWindow(:,4)))/60);

% Channel, time window
integralDiff(chanIdx,1) = abs(sum(meanSzPLI(timeWindow(:,1)))/60);
integralDiff(chanIdx,2) = abs(sum(meanSzPLI(timeWindow(:,2)))/60);
integralDiff(chanIdx,3) = abs(sum(meanSzPLI(timeWindow(:,3)))/60);
integralDiff(chanIdx,4) = abs(sum(meanSzPLI(timeWindow(:,4)))/60);


end % END FOR

% integralDiff = permute(integralDiff, [1,3,2]);

% integralDiff = reshape(integralDiff, 8,8,4);
% integralDiff = integralDiff ./ max(integralDiff(:));
%


% maxVal = 0.3211; % 2012PP05Sz1
% maxVal = 0.33;
% maxVal = 0.1937; % 2014PP04Sz4
integralDiff = integralDiff ./ maxVal;

load('E:\data\human CNS\EMD\Sz\clips\2012PP05Sz1.mat')

% integralDiff2 = integralDiff;
% integralDiff = integralDiff ./ max(integralDiff(:));

figure;
mapCol = gray(128);
border = 0;

layout = reshape(1:64,8,8);

plotData = detrend(data(10,1:sizePLI*500))./1000;

subplot(2,1,1)
hold on
box on
for kk = 1:4
    
    patchx = [timeWindowBounds(1,kk), timeWindowBounds(2,kk), timeWindowBounds(2,kk), timeWindowBounds(1,kk)];
    patchy = [2*min(plotData), 2*min(plotData), 2*max(plotData), 2*max(plotData)];
    
    pData = patch(patchx, patchy,1);
    
    set(pData, 'EdgeColor', 'none')
    set(pData, 'FaceColor', [0.5,0.5,0.5])
    set(pData, 'FaceAlpha', 0.25)
    
end % END FOR


plot(linspace(-5,5,sizePLI*500), plotData)
title('Detrended Data from Electrode 10', 'FontSize', 25)
xlabel('Time, minutes', 'FontSize', 25)
ylabel('Voltage, mV', 'FontSize', 25)
ylim([1.1*min(plotData), 1.1*max(plotData)])

a = get(gca, 'xticklabels');
b = get(gca, 'yticklabels');
set(gca, 'xticklabels', a, 'FontSize', 20)
set(gca, 'yticklabels', b, 'FontSize', 20)



for ii = 1:4
   
    subplot(2,4,ii+4)
    
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
        
        colorPatch = mapCol(floor(integralDiff(jj,ii) * 127+1),:);
        
        pData  = patch( xSubPos, ySubPos, colorPatch);
        
        
        xlim([1,size(layout,1)+1])
        ylim([1,size(layout,2)+1])
        
        set(pData, 'EdgeColor', 'k')
        axis square
        
        titleString = {'Pre-Ictal', 'Ictal', 'Early Post-Ictal', 'Late Post-Ictal'};
        tmpTitle = sprintf([titleString{ii}, '\n', num2str(timeWindowBounds(1,ii)), ' to ', num2str(timeWindowBounds(2,ii)), ' minutes']);
        title(tmpTitle, 'FontSize', 25) 
        
        set(gca, 'xtick', [1.5:8.5])
        set(gca, 'xticklabel', [1:8], 'FontSize', 20)
        
        set(gca, 'ytick', [1.5:8.5])
        set(gca, 'yticklabel', [1:8:57], 'FontSize', 20)
        
        xlabel('Channel')
        ylabel('Channel')
        box on
        
    end % END FOR num LL corners    
end

set(gcf, 'units', 'inches')
set(gcf, 'Position', [2    2   16    8])
set(gcf, 'PaperPosition', [2    2   16    8])


% EOF