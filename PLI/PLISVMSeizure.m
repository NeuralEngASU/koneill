% PLI SVM for Seizure Detection

%% Load data

%% Compute PLI

inputPLISz = 'E:\data\human CNS\EMD\Sz\ProcData\CAR';
inputPLINonSz = 'E:\data\human CNS\EMD\NonSz\ProcData\CAR';
outputPLISz = 'E:\data\PLI\EMDData\Sz';
outputPLINonSz = 'E:\data\PLI\EMDData\NonSz';

szList = dir([inputPLISz, '\*.mat']);
szFiles = szList.name;

params.winSize = 1;
params.Fs = 500;
params.chanProcess = [1:64];
params.surrFlag = 0;
params.surrNum = 0;
params.rawPhiFlag = 0;
params.biPolarFlag = 0;
params.statsFlag = 0;
params.globalFlag = 0;
params.globalChan = [0];

for ii = 1:size(szList,1)
    [~] = GenPLIChan(fullfile(inputPLISz, szList(ii).name), outputPLISz, params);
end % END FOR each file

nonszList = dir([inputPLINonSz, '\*.mat']);
nonszFiles = nonszList.name;

for jj = 1:size(nonszList,1)
    [~] = GenPLIChan(fullfile(inputPLINonSz, nonszList(jj).name), outputPLINonSz, params);
end % END FOR each file


%% Loop over all saved files
inputPLISz = 'E:\data\PLI\EMDData\Sz';
inputPLINonSz = 'E:\data\PLI\EMDData\NonSz';

rng(1);

useP = 1; % Creates virtual channels for P
useR = 1; % Creates virtual channels for R

trainRatio = 0.5; % Ratio of trials used for training

tmpData = [];

filesTrain = 3:5;%1:2;
filesTest = 3:5;

% Sz Data
for ii = filesTrain
% Load PLI
szList = dir([inputPLISz, '\*.mat']);
% szFiles = szList.name;

load(fullfile(inputPLISz, szList(1).name))
tmpData = [tmpData, [p(301:420,:)'; r(301:420,:)']];
end % Sz

trainData = tmpData;
tmpData = [];
classLabels = ones(size(tmpData,2),1);

% NonSz Data
for jj = filesTrain
% Load PLI
nonszList = dir([inputPLINonSz, '\*.mat']);
% szFiles = szList.name;

load(fullfile(inputPLINonSz, nonszList(1).name))
tmpData = [tmpData, [p(301:420,:)'; r(301:420,:)']];
end % END FOR NonSz
trainData = [trainData,tmpData];
classLabels = [classLabels; zeros(size(tmpData,2),1)];

classLabels = [ones(360,1); zeros(360,1)];

%% Randomize Columns and Rows

trainData = trainData';
rowSortTrain = [1:size(trainData,1)];
colSort = [1:size(trainData,2)];

rowSortTrain = rowSortTrain(randperm(length(rowSortTrain)));
colSort = colSort(randperm(length(colSort)));

classLabels = classLabels(rowSortTrain);
trainData = trainData(rowSortTrain,:);
trainData = trainData(:,colSort);

trainData = trainData.^2;

trainData = trainData';
%% Parameter Selection
% Automatic Cross Validation 
% Parameter selection using n-fold cross validation

stepSize = 10;
bestLog2c = 1;
bestLog2g = -1;
epsilon = 0.005;
bestcv = 0;
Ncv = 3; % Ncv-fold cross validation cross validation
deltacv = 10^6;

while abs(deltacv) > epsilon
    bestcv_prev = bestcv;
    prevStepSize = stepSize;
    stepSize = prevStepSize/2;
    log2c_list = bestLog2c-prevStepSize:stepSize/2:bestLog2c+prevStepSize;
    log2g_list = bestLog2g-prevStepSize:stepSize/2:bestLog2g+prevStepSize;
    
    numLog2c = length(log2c_list);
    numLog2g = length(log2g_list);
    cvMatrix = zeros(numLog2c,numLog2g);
    
    for i = 1:numLog2c
        log2c = log2c_list(i);
        for j = 1:numLog2g
            log2g = log2g_list(j);
            cmd = ['-q -c ', num2str(2^log2c), ' -g ', num2str(2^log2g)];
            cv = get_cv_ac(classLabels, trainData', cmd, Ncv);
            if (cv >= bestcv),
                bestcv = cv; bestc = 2^log2c; bestg = 2^log2g;
            end
        end
    end
    deltacv = bestcv - bestcv_prev;
    
end
disp(['The best parameters, yielding Accuracy=',num2str(bestcv*100),'%, are: C=',num2str(bestc),', gamma=',num2str(bestg)]);

% Best params for full time series (43.5% accuracy): bestc =  2048; bestg = 0.0028;

%% Mold Data

%% Train SVM

%% Test SVM

% cmd = ['-q -c ', num2str(2^bestc), ' -g ', num2str(2^bestg)];
cmd = ['-q -c ', num2str(bestc), ' -g ', num2str(2^bestg)];
svmModelP  = ovrtrain(classLabels,  trainData', cmd);

%% SVM Predict (LibSVM)

[label_out, acc_out, decv_out] = ovrpredictBot(classLabels, trainData', svmModelP);

% EOF