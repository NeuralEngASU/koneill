function [filePathOut] = GenPLIEpilepsy(filePath, pathOutName, params)

%% Parse Input
if ~isfield(params, 'winSize');     winSize     = 1;      else winSize     = params.winSize; end
if ~isfield(params, 'Fs');          Fs          = 500;    else Fs          = params.Fs; end
if ~isfield(params, 'chanProcess'); chanProcess = [1:32]; else chanProcess = sort(unique(params.chanProcess(:)), 'ascend'); end
if ~isfield(params, 'surrFlag');    surrFlag    = 0;      else surrFlag    = params.surrFlag; end
if ~isfield(params, 'surrNum');     surrNum     = 100;    else surrNum     = params.surrNum; end
if ~isfield(params, 'rawPhiFlag');  rawPhiFlag  = 0;      else rawPhiFlag  = params.rawPhiFlag; end
if ~isfield(params, 'biPolarFlag'); biPolarFlag = 0;      else biPolarFlag = params.biPolarFlag; end
if ~isfield(params, 'statsFlag');   statsFlag   = 0;      else statsFlag   = params.statsFlag; end
if ~isfield(params, 'globalFlag');  globalFlag  = 0;      else globalFlag  = params.globalFlag; end
if ~isfield(params, 'globalChan');  globalChan  = [1:64]; else globalChan  = sort(unique(params.globalChan(:)), 'ascend'); end

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

if globalFlag
    filePathOut = [filePathOut, '_Global'];
    if sum(diff(globalChan)==1) == size(globalChan(:), 1)-1
        filePathOut = [filePathOut, sprintf('%d-%d', min(globalChan), max(globalChan))];
    else
        filePathOut = [filePathOut, sprintf('%d~%d', min(globalChan), max(globalChan))];
    end % END IF isConsecutive
    
end % END IF globalFlag

if biPolarFlag
    filePathOut = [filePathOut, '_BiPolar'];
end % END IF biPolarFlag

if rawPhiFlag
    filePathOut = [filePathOut, '_DeltaPhi'];
end % END IF statsFlag

if statsFlag
    filePathOut = [filePathOut, '_Stats'];
end % END IF statsFlag

% filePathOut = fullfile(pathOutName, [fileName, '_PLI_winSize', num2str(winSize), '_rawPhi.mat']);
% filePathOut = fullfile(pathOutName, [fileName, '_PLI_winSize', num2str(winSize), '.mat']);
filePathOut = [filePathOut, '_PLI_winSize', num2str(winSize), '.mat'];

%% Give User Feedback
fprintf('******************************************************************\n')
fprintf('Generate Phase Lag Index Data\n')
fprintf('Loading Data: %s\n', filePath)
fprintf('Data Size: %f MB\n',sizeFile)

%% Parse variables
load(filePath);
% Header = openNSx(filePath);
% Header = Header.MetaTags;
% Header = header;

% Detrend Data

if size(data,1) > size(data,2)
    data = detrend(data);
else
    data = detrend(data');
end % END IF

% Find Instantaneous Phase
data = atan2(imag(hilbert(data)),data);

if isempty(chanProcess)
    chanProcess = [1:size(data,2)]';
end %END IF
    


%%

% Find the number of channels
% numChans = size(chanProcess,1);
numChans = size(data,2);
% chanProcess = [1:numChans]';
params.chanProcess = chanProcess;

% Set up channel parings
chanPairNums = nchoosek(sort(unique(1:numChans),'ascend'),2);
pairNum      = size(chanPairNums,1);

% Calculate Number of Windows
% dataLen = Header.DataPoints;
dataLen = size(data,1);
winNum  = floor(dataLen / (winSize * Fs));
%%

% Seperate Data into Windows
% Windowed data x winNum x channel
rawWin = reshape(permute(data(1:floor(winNum*winSize*Fs),:), [1,3,2]), floor(winSize * Fs),winNum,size(chanProcess,1));

deltaPhi = 0;
smp = 0;

% initialize variables
if surrFlag
    % Set up PLI and R variables
    pli = zeros(winNum, pairNum, surrNum+1);
    r   = zeros(winNum, pairNum, surrNum+1);
    
    % for surrogate data: random phase shift applied to each channel
    smp = [zeros(pairNum,1) round(Fs)/2 + round(round(Fs)*rand(pairNum,surrNum))]; % 500-1000 samples (0.5 - 1.0 seconds)
    
    if rawPhiFlag
        deltaPhi = zeros(winNum, pairNum, floor(winSize*Fs), surrNum+1);
    end % END IF rawPhiFlag
    
else
    pli = zeros(winNum, pairNum, 1);
    
%     % Compute instantaneous phi for all channels
    for phiChan = 1:size(rawWin,3)
        phi(:,:,phiChan) = atan2(imag(hilbert(rawWin(:,:,phiChan))),rawWin(:,:,phiChan));
    end
    % stats
    r   = zeros(winNum, pairNum, 1);
    
    if statsFlag
        binEdge = [-pi:pi/100:pi];
        
        circMean = zeros(winNum, pairNum, 3); % circular mean, upper bound, lower bound
        circMed  = zeros(winNum, pairNum, 1); % Cirular median
        circVar  = zeros(winNum, pairNum, 1); % circular variance
        vMParams = zeros(winNum, pairNum, 2); % thetahat, kappa
        vMScale  = zeros(winNum, pairNum, 1); % Scaling params for vM (transfers the units to [counts])
        RMSE     = zeros(winNum, pairNum, 2); % Mean squared Error
        vMCorr   = zeros(winNum, pairNum, 1); % Circular and Angular Correlation coefficiant between von Mises and hist(delaPhi)
        vMR2     = zeros(winNum, pairNum, 1); % R-squared value
        circStd  = zeros(winNum, pairNum, 2); % angular and circular std
        circSkew = zeros(winNum, pairNum, 2); % Pewsey and Fischer Skewness
        circKurt = zeros(winNum, pairNum, 2); % Pewsey and Fischer kurtosis
    end % END IF statsFlag
    
    if rawPhiFlag
        deltaPhi = zeros(winNum, pairNum, floor(winSize*Fs), 1);
    end % END IF rawPhiFlag
    
end % END IF surrFlag

if exist('header', 'var')
    Header = header;
    Header.FilePathIn = pathInName;
    Header.FileNameIn = fileName;
    Header.FileOutLoc = filePathOut;
else
    Header.Fs = Fs;
    Header.FilePathIn = pathInName;
    Header.FileNameIn = fileName;
    Header.FileOutLoc = filePathOut;
end % END IF exist(header)

Header.params = params;
Header.params.winNum = winNum;
Header.chanPairNums = chanPairNums;

%% Set up parrallel workers

if surrFlag
    parpool(8);
end % END IF surrFlag

%% Proccess Data
fprintf('*****\nProcessing Data\n')
timeWatch = tic;

% Create save file
save(filePathOut,'chanPairNums','Header', 'params', 'phi','-v7.3');

% For each channel pair
for cp = 1:pairNum
    
    % pull out data for this iteration
    % windowed Data x winNum x channel
    raw1 = double(rawWin(:,:,chanPairNums(cp,1)));
    raw2 = double(rawWin(:,:,chanPairNums(cp,2)));
    
    if surrFlag
         % Create temporary PLI and R data
         tmpp = nan(winNum, surrNum+1);
         tmpr = nan(winNum, surrNum+1);
         tmpsmp = smp(cp,:);
        
        parfor ss = 1:surrNum+1
            
            % circshift the second channel (destroy correlations)
            raw2s = circshift(raw2,[0, tmpsmp(ss), 0]);
            
            % calculate PLI and R
            if rawPhiFlag
%                 tmpDeltaPhi = nan(winNum, (winSize*Fs), surrNum+1);
                [tmpp(:,ss),tmpr(:,ss), tmpDeltaPhi(:,:,ss)]=pli(raw1, raw2s);
            else
                [tmpp(:,ss),tmpr(:,ss),~]=pli(raw1, raw2s);
            end % END IF rawPhiFlag
        
        end
    else % No surrogates
        
        %%
        tmpp = nan(winNum, 1);
        tmpr = nan(winNum, 1);
        
        if rawPhiFlag
            tmpDeltaPhi = nan(winNum, 1,(winSize*Fs));
            [tmpp(:,1),tmpr(:,1), tmpDeltaPhi(:,1,:)] = wpliECoG(raw1, raw2);
        else
            [tmpp(:,1),tmpr(:,1)] = wpliEpilepsy(raw1, raw2);
        end % END IF rawPhiFlag
        
        if statsFlag
            [tmpCircMean, tmpCircMed, tmpCircVar, tmpvMParams,...
             tmpvMScale, tmpRMSE , tmpvMCorr, tmpvMR2, ...
             tmpCircStd,tmpCircSkew, tmpCircKurt] = CircStats(squeeze(tmpDeltaPhi)');
        end % END IF statsFlag   
     %%
    end % END IF surrFlag
    
    p(:,cp,:) = tmpp;
    r(:,cp,:) = tmpr;
%     eval(['P', num2str(cp), ' = tmpp;']);
%     eval(['R', num2str(cp), ' = tmpp;']);
    
    if statsFlag
        circMean(:,cp,:) = tmpCircMean;
        circMed(:,cp,:)  = tmpCircMed;
        circVar(:,cp,:)  = tmpCircVar;
        vMParams(:,cp,:) = tmpvMParams;
        vMScale(:,cp,:)  = tmpvMScale;
        RMSE(:,cp,:)     = tmpRMSE;
        vMCorr(:,cp,:)   = tmpvMCorr;
        vMR2(:,cp,:)     = tmpvMR2;
        circStd(:,cp,:)  = tmpCircStd;
        circSkew(:,cp,:) = tmpCircSkew;
        circKurt(:,cp,:) = tmpCircKurt;
    end % END IF statsFlag
    
    
    if rawPhiFlag
        deltaPhi(:,cp,:,:) = tmpDeltaPhi;
    end % END IF rawPhiFlag
    
    timeSpent = toc(timeWatch);
    
    if surrFlag
        fprintf('Channel Pair %d/%d\n', cp, pairNum)
    else
        if ~mod(cp,1)
            clc;
            fprintf('Channel Pair: %d/%d\n', cp, pairNum)
            fprintf('Time Spent: %f\n', timeSpent)
            fprintf('Time Spent Per Pair: %f\n', timeSpent/cp)
            fprintf('Time Left: %f\n', (timeSpent/cp)*(pairNum-cp))
        end % END IF
    end % END IF surrFlag
    
% save(filePathOut,['P', num2str(cp)], ['R', num2str(cp)], '-append');
    
end % END FOR

fprintf('Saving Data to: %s\n', filePathOut)
% save results

% if globalFlag
%     chanPairNums = chanPairNums-1;
% end % END IF globalFlag

% save(filePathOut,'p','r','chanPairNums','smp', 'Header', 'deltaPhi', '-v7.3');
save(filePathOut,'p','r', '-append');
% if globalFlag
%     save(filePathOut, 'globalChan', '-append')
% end % END IF globalFlag
% 
% if statsFlag
%     [~, binCenter] = hist(ones(1,size(binEdge,1)),binEdge);
%     
%     save(filePathOut, 'circMean', 'circMed', 'circVar','vMParams',...
%                       'vMScale','RMSE','vMCorr','vMR2','circStd',...
%                       'circSkew','circKurt', 'binEdge', 'binCenter', '-append')
% end % END IF statsFlag
end % END FUNCTION
% EOF