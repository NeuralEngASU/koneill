function [ Header ] = GenCohereEpilepsy(filePath, pathOutName, params)

%% Parse Input  
if ~isfield(params, 'winSize');     winSize            = 1;        else winSize            = params.winSize; end
if ~isfield(params, 'Fs');          Fs                 = 500;      else Fs                 = params.Fs; end
if ~isfield(params, 'chanProcess'); chanProcess        = [1:32];   else chanProcess        = sort(unique(params.chanProcess(:)), 'ascend'); end
if ~isfield(params, 'sizeWindow');  cParams.sizeWindow = Fs;       else cParams.sizeWindow = params.sizeWindow; end
if ~isfield(params, 'numOverlap');  cParams.numOverlap = 1;        else cParams.numOverlap = params.numOverlap; end
if ~isfield(params, 'freqBand');    cParams.freqBand   = 1:(Fs/2); else cParams.freqBand   = params.freqBand; end
if ~isfield(params, 'numWindow');   numWindow          = 3;        else numWindow          = params.numWindow; end

cParams.Fs = Fs;

% Parse file strings
pathInExpr = '(.+)?\\';
pathInName = regexp(filePath, pathInExpr, 'Tokens');
pathInName = pathInName{1}{1};

nameExpr = '.+\\(.+\..+$)';
fileName = regexp(filePath, nameExpr, 'Tokens');
fileName = fileName{1}{1};

sizeFile = dir(filePath);
sizeFile = round(sizeFile.bytes / 1024^2, 1);

filePathOut = fullfile(pathOutName, fileName(1:end-4));

filePathOut = [filePathOut, '_Cohere_winSize', num2str(winSize), '_numWindow', num2str(numWindow), '.mat'];

%% Give User Feedback
fprintf('******************************************************************\n')
fprintf('Generate Coherence Index Data\n')
fprintf('Loading Data: %s\n', filePath)
fprintf('Data Size: %f MB\n',sizeFile)

%% Parse variables

load(filePath);

if size(data,1) > size(data,2)
    data = detrend(data);
else
    data = detrend(data');
end % END IF

if isempty(chanProcess)
    chanProcess = [1:size(data,2)]';
end %END IF

%%
numChans = size(data,2);
params.chanProcess = chanProcess;

% Set up channel parings
chanPairNums = nchoosek(sort(unique(1:numChans),'ascend'),2);
pairNum      = size(chanPairNums,1);

% Calculate Number of Windows
% dataLen = Header.DataPoints;
dataLen = size(data,1);
winNum  = floor(dataLen / (winSize * Fs));

%% Reshape data  

numWindow = 3;
boundErr = ceil(numWindow/2);
distWindow = [boundErr-numWindow, -(boundErr-numWindow)];

rawWin = zeros(cParams.sizeWindow*numWindow, winNum, numChans);

for win = 1:winNum
    if win < boundErr 
        dataIdx = [1:(cParams.sizeWindow*numWindow)];
        rawWin(:,win,:) = permute(data(dataIdx,:), [1,3,2]);
    elseif win > winNum - boundErr
        dataIdx = [1:(cParams.sizeWindow*numWindow)] + (winNum-boundErr-1)*cParams.sizeWindow;
        rawWin(:,win,:) = permute(data(dataIdx,:), [1,3,2]);
    else
        dataIdx = [1:(cParams.sizeWindow*numWindow)] + (win-boundErr)*cParams.sizeWindow;
        rawWin(:,win,:) = permute(data(dataIdx,:), [1,3,2]);
    end
end % END FOR each window

%% Proccess Data
fprintf('*****\nProcessing Data\n')
timeWatch = tic;

c = zeros(winNum, pairNum);

% Create save file
save(filePathOut,'chanPairNums','Header', 'params','-v7.3');
parpool(16)
parfor cp = 1:pairNum
    
    % pull out data for this iteration
    % windowed Data x winNum x channel
    raw1 = double(rawWin(:,:,chanPairNums(cp,1)));
    raw2 = double(rawWin(:,:,chanPairNums(cp,2)));
    
     c(:,cp) = cohereEpilepsy(raw1, raw2, cParams);
    
%     c(:,cp) = tmpc;
    
    timeSpent = toc(timeWatch);
    
    if ~mod(cp,1)
        clc;
        fprintf('Channel Pair: %d/%d\n', cp, pairNum)
        fprintf('Time Spent: %f\n', timeSpent)
        fprintf('Time Spent Per Pair: %f\n', timeSpent/cp)
        fprintf('Time Left: %f\n', (timeSpent/cp)*(pairNum-cp))
    end % END IF
    
end % END FOR each channel pair

poolobj = gcp('nocreate');
delete(poolobj);

fprintf('Saving Data to: %s\n', filePathOut)

% save results
save(filePathOut,'c','-append');

end % END FUNCTION

% EOF