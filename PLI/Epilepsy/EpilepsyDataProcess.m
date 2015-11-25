%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Epilepsy Data Process
%   Lab: Neural Engineering Laboratory
%   PI: Bradley Greger
%   Author: Kevin O'Neill
%   Date: 2015.11.12
%
%   Desc: Beginning-to-end processing of epilepsy data.
%       - Filter (notch filter and band-pass (optional)
%       - Generate WPLI
%       - Load GridDef (must be manually coded in GridDef.m)
%       - Compute Neighbor-Average
%       - Plot the time-evolution for each channel (positons based on
%       GridDef.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Chosen Files
tic
filePath = 'E:\data\human CNS\EMD\Sz\clips\';
fileName = '2014PP04Sz5';

%% Freq Process
clear Header
Header.rawFileName = fileName;
Header.rawFilePath = filePath;
Header.Fs = 500;


params.targetDir = 'D:\PLI\SeizureDetection\PreProcessedData\Sz\';
params.combFilter = 1;
params.combFilterBand = []; % Default to 60:60:(Fs/2)
params.bandPassFilter = 1;
params.passBands = [0,30; 30,80; 80,250];

[fileList] = EpilepsyFreqProcess(filePath, [fileName, '.mat'], params, Header);

%% Generate PLI

filePath = 'D:\PLI\SeizureDetection\PreProcessedData\Sz\';
% fileName = fileList;

% test = dir('D:\PLI\SeizureDetection\PreProcessedData\Sz');
% test.name


for ii = 1:length(fileList)
    params.targetDir = 'D:\PLI\SeizureDetection\PLIProcessedData\Sz\';
    
    params.winSize = 1;
    params.Fs = 500;
    params.chanProcess = [];
    params.surrFlag = 0;
    params.surrNum = 0;
    params.rawPhiFlag = 0;
    params.biPolarFlag = 0;
    params.statsFlag = 0;
    params.globalFlag = 0;
    params.globalChan = [];
    
    [~] = GenPLIEpilepsy(fullfile(filePath, fileList{ii}), params.targetDir, params);

end % END FOR

toc
%% Generate Coherence (mscohere, cpsd)

filePath = 'D:\PLI\SeizureDetection\PreProcessedData\Sz\';

files = dir(fullfile(filePath, [fileName, '*.mat']));

fileList = {files.name};

expr = '([0-9])+_([0-9])+\.mat';

poolobj = gcp('nocreate');
delete(poolobj);

for ii = 1:length(fileList)
    
    
    test = regexp(fileList{ii}, expr, 'tokens');
    
    if isempty(test)
        params.freqBand = 1:(Fs/2);
    else
        params.freqBand = str2double(test{1}{1}) : str2double(test{1}{2});
    end % END IF freqBand
        
    params.targetDir = 'D:\PLI\SeizureDetection\CohereProcessedData\Sz';
    
    params.winSize = 1;
    params.Fs = 500;
    params.chanProcess = [1:32];
    params.sizeWindow = Fs;
    params.numOverlap = 1;
%     params.freqBand = 1:(Fs/2);
    params.numWindow = 3;
    
    [~] = GenCohereEpilepsy(fullfile(filePath, fileList{ii}), params.targetDir, params);

end % END FOR



%% Load GridDef

%% Neighbor Computation

%% Plot Time-Series
% fileName = '2014PP04NonSz4';
% fileLoad = fileName(1:end-4);
filePathProc   = 'D:\PLI\SeizureDetection\PreProcessedData\Sz\';
filePathPLI    = 'D:\PLI\SeizureDetection\PLIProcessedData\Sz\';
filePathCohere = 'D:\PLI\SeizureDetection\CohereProcessedData\Sz\';

files = dir(fullfile(filePathProc,   [fileName, '*.mat']));
fileListProc = {files.name};
files = dir(fullfile(filePathPLI,    [fileName, '*.mat']));
fileListPLI = {files.name};
files = dir(fullfile(filePathCohere, [fileName, '*.mat']));
fileListCohere = {files.name};


load([filePathPLI, fileListPLI{end}])

chanIdx = 10;

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

idx = zeros(size(chanPairNums,1),1);
% Find idicies
for jj = 1:size(desiredChanPairs,1)
    idx = idx | ((chanPairNums(:,1) == desiredChanPairs(jj,1)) & (chanPairNums(:,2) == desiredChanPairs(jj,2)));
end

% idx = find(idx==1);
% idx = false(size(chanPairNums,1),1);
% idx = 7140;
% chanTest = [82, 112];
% idx = idx | ((chanPairNums(:,1) == chanTest(1)) & (chanPairNums(:,2) == chanTest(2)));

idx = find(idx==1);

t = linspace(-5,5,size(p,1));
tRaw = linspace(-5,5, size(p,1)*Fs);

figure

titleStr = {['Comb ', fileName], '0:30 Hz', '30:80 Hz', '80:250 Hz'};

for ii = 1:4

    plotIdx = [1,4,7,10];
    fileIdx = [4,1,2,3];
    
    load([filePathProc, fileListProc{fileIdx(ii)}]);
    load([filePathPLI, fileListPLI{fileIdx(ii)}]);
    load([filePathCohere, fileListCohere{fileIdx(ii)}]);
    
    meanPLI = mean(p(:,idx),2);
    meanCohere = mean(c(:,idx),2);
    
    stdErrPLI = 2*(std(p(:,idx),0,2)/sqrt(length(idx)));
    stdErrCohere = 2*(std(c(:,idx),0,2)/sqrt(length(idx)));
    
    
    ax1(ii) = subplot(4,3,plotIdx(ii));
    plot(tRaw, data(1:(size(p,1)*Fs),chanIdx))
    ylabel('Voltage')
    title([titleStr{ii}, ' Processed Data'])
    
    ax2(ii) = subplot(4,3,plotIdx(ii)+1);
    pData = patch([t, fliplr(t)],[meanPLI+stdErrPLI; flipud(meanPLI-stdErrPLI)], 'k');
    hold on
    plot(t,meanPLI)
    hold off
    title([titleStr{ii}, ' PLI Data'])
    ylabel('PLI')
    set(pData, 'FaceAlpha', 0.25)
    set(pData, 'EdgeColor', 'none')
    
    ax3(ii) = subplot(4,3,plotIdx(ii)+2);
    pData = patch([t, fliplr(t)],[meanCohere+stdErrCohere; flipud(meanCohere-stdErrCohere)], 'k');
    hold on
    plot(t,meanCohere)
    hold off
    title([titleStr{ii}, ' Coherence Data'])
    ylabel('Coherence')
    set(pData, 'FaceAlpha', 0.25)
    set(pData, 'EdgeColor', 'none')
    
end % END FOR

subplot(4,3,10)
xlabel('Time, minutes')

subplot(4,3,11)
xlabel('Time, minutes')

subplot(4,3,12)
xlabel('Time, minutes')

linkaxes([ax1,ax2,ax3],'x')
toc
%%
% 
% Fs = 500;  % Sampling Frequency
% 
% Fpass1 = 54;          % First Passband Frequency
% Fstop1 = 59;          % First Stopband Frequency
% Fstop2 = 61;          % Second Stopband Frequency
% Fpass2 = 66;          % Second Passband Frequency
% Apass1 = 1;           % First Passband Ripple (dB)
% Astop  = 60;          % Stopband Attenuation (dB)
% Apass2 = 1;           % Second Passband Ripple (dB)
% match  = 'stopband';  % Band to match exactly
% 
% % Construct an FDESIGN object and call its BUTTER method.
% h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
%                       Apass2, Fs);
% Hd = design(h, 'butter', 'MatchExactly', match);
% 
% for ii = 1:size(data,2)
%    dataFilt(:,ii) = filtfilt(Hd.sosMatrix, Hd.ScaleValues, data(:,ii));
%    disp(ii)
% end





% 
% % EOF