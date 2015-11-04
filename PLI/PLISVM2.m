%% PLI SVM Time Series and Subsets
% Subset 1: Facemotor grid coupling
% Subset 2: Wernekies Grid Coupling
% Subest 3: Inter-Grid Coupling
% Subest 4-?: Frequency band Coupling for above subsets

%%
load('E:\data\PLI\delta\PLIOutput\Delta_ProcessedTrialData_PLI_winSize0.1.mat')

trainRatio = 0.5; % The ratio of trials to train from

uniqueTrials = sort(unique(Header.class), 'ascend');

rng(1); % For reproducibility

count = size(Header.class,2);
maxTrials = count; % The maximum number of trials to use. This ensures no class is over represented

for ii = 1:size(uniqueTrials,2)
    count = sum(Header.class == uniqueTrials(ii));
    if maxTrials > count
        maxTrials = count;
    end % END IF find max trials
end % END FOR find max trials


for ii = 1:size(uniqueTrials,2)
    
    logicalMat = logical(zeros(1, size(Header.class,2)));
%     numTrials = sum(Header.class == uniqueTrials(ii));
    trialIdx = find(Header.class == uniqueTrials(ii));
    
    randIdx = randperm(size(trialIdx,2));
    
    floorVal = floor(maxTrials*trainRatio);
    ceilVal  = ceil(maxTrials*trainRatio);
    
    if floorVal == ceilVal
        ceilVal = ceilVal + 1;
    end % END IF
    
    trainIdx = trialIdx(randIdx(1:floorVal));
    testIdx = trialIdx(randIdx(ceilVal:ceilVal+floor(maxTrials/2)));
    
    trainTrial{ii} = logicalMat;
    testTrial{ii}  = logicalMat;
    
    trainTrial{ii}(trainIdx) = 1;
    testTrial{ii}(testIdx) = 1;
        
end % END FOR uniqueTrials

%% Select and mold data

for jj = 1:size(uniqueTrials,2)
    
    tmpDataP = squeeze(p(11:20,:,:,trainTrial{jj}));
    tmpDataP = permute(tmpDataP, [2,1,3]);
    tmpDataP = reshape(tmpDataP, size(tmpDataP,1), 10*size(tmpDataP,3));
    
    trainDataP{jj}  = tmpDataP;
    
    tmpDataP = squeeze(p(1:10,:,:,trainTrial{jj}));
    tmpDataP = permute(tmpDataP, [2,1,3]);
    tmpDataP = reshape(tmpDataP, size(tmpDataP,1), 10*size(tmpDataP,3));

    trainZeroP{jj} = tmpDataP;


    tmpDataP = squeeze(p(11:20,:,:,testTrial{jj}));
    tmpDataP = permute(tmpDataP, [2,1,3]);
    tmpDataP = reshape(tmpDataP, size(tmpDataP,1), 10*size(tmpDataP,3));

    testDataP{jj}  = tmpDataP;

    tmpDataP = squeeze(p(1:10,:,:,testTrial{jj}));
    tmpDataP = permute(tmpDataP, [2,1,3]);
    tmpDataP = reshape(tmpDataP, size(tmpDataP,1), 10*size(tmpDataP,3));
    
    testZeroP{jj} = tmpDataP;
    
end % END FOR uniqueTrials

sizeTrainClass = 0;

trainMapP  = [];
classMapTrainP   = [];

testMapP  = [];
classMapTestP   = [];

trainZeroMapP = [];
testZeroMapP = [];

for jj = 1:size(uniqueTrials,2)
    
    trainMapP  = [trainMapP,  trainDataP{jj} ];

    classMapTrainP   = [classMapTrainP;  ones(size(trainDataP{jj},2),1) * jj ];
 
    testMapP  = [testMapP,  testDataP{jj} ];
   
    classMapTestP   = [classMapTestP;  ones(size(testDataP{jj},2),1) * jj ];
   
    trainZeroMapP = [trainZeroMapP, trainZeroP{jj}]; 
    testZeroMapP = [testZeroMapP, testZeroP{jj}]; 

end % END FOR uniqueTrials

% Insert 0 class
sizeTrainZero = size(trainZeroMapP,2);
sizeTestZero  = size(testZeroMapP,2);

sizeTrainClass  = size(trainDataP{1},2);
sizeTestClass   = size(testDataP{1},2);

availTrainZero = randperm(sizeTrainZero);
availTestZero  = randperm(sizeTestZero);

trainMapP = [trainMapP, trainZeroMapP(:,availTrainZero(1:sizeTrainClass))];
testMapP  = [testMapP, testZeroMapP(:,availTestZero(1:sizeTestClass))];

classMapTrainP = [classMapTrainP; ones(sizeTrainClass,1)*11];
classMapTestP = [classMapTestP; ones(sizeTestClass,1)*11];

% Randomly sort the rows and columns
testMapP = testMapP';
trainMapP = trainMapP';

rowSortTrain = [1:size(trainMapP,1)];
rowSortTest = [1:size(testMapP,1)];
colSort = [1:size(trainMapP,2)];

rowSortTrain = rowSortTrain(randperm(length(rowSortTrain)));
rowSortTest = rowSortTest(randperm(length(rowSortTest)));
colSort = colSort(randperm(length(colSort)));

classMapTrainP = classMapTrainP(rowSortTrain);
trainMapP = trainMapP(rowSortTrain,:);
trainMapP = trainMapP(:,colSort);

classMapTestP = classMapTestP(rowSortTest);
testMapP = testMapP(rowSortTest,:);
testMapP = testMapP(:,colSort);

trainMapP  = trainMapP.^2;

testMapP  = testMapP.^2;

testMapP = testMapP';
trainMapP = trainMapP';

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
            cv = get_cv_ac(classMapTrainP, trainMapP', cmd, Ncv);
            if (cv >= bestcv),
                bestcv = cv; bestc = 2^log2c; bestg = 2^log2g;
            end
        end
    end
    deltacv = bestcv - bestcv_prev;
    
end
disp(['The best parameters, yielding Accuracy=',num2str(bestcv*100),'%, are: C=',num2str(bestc),', gamma=',num2str(bestg)]);

% Best params for full time series (43.5% accuracy): bestc = 0.0020; bestg = 0.0156;
% Best params for full Mean/Diff/Max (41.6% accuracy): C = 0.35355; gamma = 0.0027621;

%% SVM Train and Test (MATLAB)

t = templateSVM('Standardize', 1, 'KernelFunction', 'gaussian', 'BoxConstraint',bestc, 'KernelScale', bestg);
Mdl = fitcecoc(trainMapP', classMapTrainP, 'Learners', t, 'Coding', 'onevsone');


% SVM Predict (MATLAB)
predicted = predict(Mdl, testMapP');
acc = sum(predicted == classMapTestP)/numel(predicted)
