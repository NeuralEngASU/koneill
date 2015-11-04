% PLI SVM
% Use SVM to classify PLI

%% Seperate the training and test group
load('E:\data\PLI\delta\PLIOutput\Delta_ProcessedTrialData_PLI_winSize0.1.mat')

trainRatio = 0.5; % The ratio of trials to train from

uniqueTrials = sort(unique(Header.class), 'ascend');

rng(1); % For reproducibility

for ii = 1:size(uniqueTrials,2)
    
    logicalMat = false(1, size(Header.class,2));
    numTrials = sum(Header.class == uniqueTrials(ii));
    trialIdx = find(Header.class == uniqueTrials(ii));
    
    randIdx = randperm(size(trialIdx,2));
    
    floorVal = floor(numTrials*trainRatio);
    ceilVal  = ceil(numTrials*trainRatio);
    
    if floorVal == ceilVal
        ceilVal = ceilVal + 1;
    end % END IF
    
    trainIdx = trialIdx(randIdx(1:floorVal));
    testIdx = trialIdx(randIdx(ceilVal:end));
    
    trainTrial{ii} = logicalMat;
    testTrial{ii}  = logicalMat;
    
    trainTrial{ii}(trainIdx) = 1;
    testTrial{ii}(testIdx) = 1;
        
end % END FOR uniqueTrials

%% Select and mold data

for jj = 1:size(uniqueTrials,2)
    
%     tmpDataP = squeeze(p(:,:,:,trainTrial{jj})).^2;
%     tmpDataP = tmpDataP(11:19,:,:);
%     meanDataP = mean(tmpDataP);
%     diffDataP = mean(diff(tmpDataP));
%     maxDataP = max(tmpDataP);
%     tmpDataP = [meanDataP;diffDataP;maxDataP];
%     tmpDataP = reshape(tmpDataP, size(tmpDataP,1)*size(tmpDataP,2), 1,size(tmpDataP,3));
%     tmpDataP = reshape(tmpDataP, size(tmpDataP,1), size(tmpDataP,3))';

    tmpDataP = squeeze(p(:,:,:,trainTrial{jj}));
    tmpDataP = reshape(tmpDataP, 25*size(tmpDataP,2), 1,size(tmpDataP,3));
    tmpDataP = reshape(tmpDataP, size(tmpDataP,1), size(tmpDataP,3))';

%     tmpDataR = squeeze(r(:,:,:,trainTrial{jj}));
%     tmpDataR = reshape(tmpDataR, 25*size(tmpDataR,2), 1,size(tmpDataR,3));
%     tmpDataR = reshape(tmpDataR, size(tmpDataR,1), size(tmpDataR,3))';
    
    trainDataP{jj}  = tmpDataP;
%     trainDataR{jj}  = tmpDataR;
%     trainDataPR{jj} = [tmpDataP, tmpDataR];
% 
%     tmpDataP = squeeze(p(:,:,:,testTrial{jj}));
%     tmpDataP = reshape(tmpDataP, 25*size(tmpDataP,2), 1,size(tmpDataP,3));
%     tmpDataP = reshape(tmpDataP, size(tmpDataP,1), size(tmpDataP,3))';
    
%     tmpDataP = squeeze(p(:,:,:,testTrial{jj})).^2;
%     tmpDataP = tmpDataP(11:19,:,:);
%     meanDataP = mean(tmpDataP);
%     diffDataP = mean(diff(tmpDataP));
%     maxDataP = max(tmpDataP);
%     tmpDataP = [meanDataP;diffDataP;maxDataP];
%     tmpDataP = reshape(tmpDataP, size(tmpDataP,1)*size(tmpDataP,2), 1,size(tmpDataP,3));
%     tmpDataP = reshape(tmpDataP, size(tmpDataP,1), size(tmpDataP,3))';

    tmpDataP = squeeze(p(:,:,:,testTrial{jj}));
    tmpDataP = reshape(tmpDataP, 25*size(tmpDataP,2), 1,size(tmpDataP,3));
    tmpDataP = reshape(tmpDataP, size(tmpDataP,1), size(tmpDataP,3))';

%     tmpDataR = squeeze(r(:,:,:,testTrial{jj}));
%     tmpDataR = reshape(tmpDataR, 25*size(tmpDataR,2), 1,size(tmpDataR,3));
%     tmpDataR = reshape(tmpDataR, size(tmpDataR,1), size(tmpDataR,3))';
% 
%     
    testDataP{jj}  = tmpDataP;
%     testDataR{jj}  = tmpDataR;
%     testDataPR{jj} = [tmpDataP, tmpDataR];
    
end % END FOR uniqueTrials

trainMapP  = [];
trainMapR  = [];
trainMapPR = [];

classMapTrainP   = [];
classMapTrainR   = [];
classMapTrainPR  = [];

testMapP  = [];
testMapR  = [];
testMapPR = [];

classMapTestP   = [];
classMapTestR   = [];
classMapTestPR  = [];

for jj = 1:size(uniqueTrials,2)
    
    trainMapP  = [trainMapP;  trainDataP{jj} ];
%     trainMapR  = [trainMapR;  trainDataR{jj} ];
%     trainMapPR = [trainMapPR; trainDataPR{jj}];
    
    classMapTrainP   = [classMapTrainP;  ones(size(trainDataP{jj},1),1) * jj ];
%     classMapTrainR   = [classMapTrainR;  ones(size(trainDataR{jj},1),1) * jj ];
%     classMapTrainPR  = [classMapTrainPR; ones(size(trainDataPR{jj},1),1) * jj];
    
    testMapP  = [testMapP;  testDataP{jj} ];
%     testMapR  = [testMapR;  testDataR{jj} ];
%     testMapPR = [testMapPR; testDataPR{jj}];
    
    classMapTestP   = [classMapTestP;  ones(size(testDataP{jj},1),1) * jj ];
%     classMapTestR   = [classMapTestR;  ones(size(testDataR{jj},1),1) * jj ];
%     classMapTestPR  = [classMapTestPR; ones(size(testDataPR{jj},1),1) * jj];
    
    
end % END FOR uniqueTrials


% Randomly sort the rows and columns
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
% trainMapR  = trainMapR.^2;
% trainMapPR = trainMapPR.^2;

testMapP  = testMapP.^2;
% testMapR  = testMapR.^2;
% testMapPR = testMapPR.^2;

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
            cv = get_cv_ac(classMapTrainP, trainMapP, cmd, Ncv);
            if (cv >= bestcv),
                bestcv = cv; bestc = 2^log2c; bestg = 2^log2g;
            end
        end
    end
    deltacv = bestcv - bestcv_prev;
    
end
disp(['The best parameters, yielding Accuracy=',num2str(bestcv*100),'%, are: C=',num2str(bestc),', gamma=',num2str(bestg)]);

% Best params for full time series (43.5% accuracy): bestc = 64; bestg = 0.00048828;
% Best params for full Mean/Diff/Max (41.6% accuracy): C = 0.35355; gamma = 0.0027621;


%% SVM Train (LibSVM)
    
% svmModelP  = svmtrain(classMapTrainP,  trainMapP);
% svmModelR  = svmtrain(classMapTrainR,  trainMapR);
% svmModelPR = svmtrain(classMapTrainPR, trainMapPR);
cmd = ['-q -c ', num2str(2^bestc), ' -g ', num2str(2^bestg)];
svmModelP  = ovrtrain(classMapTrainP,  trainMapP', cmd);
% svmModelR  = ovrtrain(classMapTrainR,  trainMapR, cmd);
% svmModelPR = ovrtrain(classMapTrainPR, trainMapPR, cmd);

%% SVM Predict (LibSVM)

% [predictClassP, accP, probP]  = svmpredict(classMapTestP,  testMapP, svmModelP);
% [predictClassP]  = svmpredict(classMapTestP,  testMapP, svmModelP);
% [predictClassR]  = svmpredict(classMapTestR,  testMapR, svmModelR);
% [predictClassPR] = svmpredict(classMapTestPR,  testMapPR, svmModelPR);

[label_out, acc_out, decv_out] = ovrpredictBot(classMapTestP, testMapP', svmModelP);

%% SVM Train (MATLAB)

t = templateSVM('Standardize', 1, 'KernelFunction', 'gaussian', 'BoxConstraint',bestc, 'KernelScale', bestg);
Mdl = fitcecoc(trainMapP', classMapTrainP, 'Learners', t, 'Coding', 'onevsone');


% SVM Predict (MATLAB)
predicted = predict(Mdl, testMapP');
acc = sum(predicted == classMapTestP)/numel(predicted)

% EOF