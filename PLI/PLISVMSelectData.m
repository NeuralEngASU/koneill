function  [trainData, trainLabels, testData, testLabels] = PLISVMSelectData(dataParams, tryRNG)

load(fullfile(dataParams.pathName, dataParams.fileName))

trainRatio = dataParams.trainRatio; % The ratio of trials to train from

uniqueTrials = sort(unique(Header.class), 'ascend');

rng(tryRNG); % For reproducibility

count = size(Header.class,2);
maxTrials = count; % The maximum number of trials to use. This ensures no class is over represented

for ii = 1:size(dataParams.classes,2)
    count(ii) = sum(Header.class == dataParams.classes(ii));
end % END FOR find max trials for each class

% Extract data for the specified classes

for ii = 1:size(dataParams.classes,2)
    
    logicalMat = false(1, size(Header.class,2));
%     numTrials = sum(Header.class == uniqueTrials(ii));
    trialIdx = find(Header.class == dataParams.classes(ii));
    
    randIdx = randperm(size(trialIdx,2));
    
    floorVal = floor(count(ii)*trainRatio);
    ceilVal  = ceil(count(ii)*trainRatio);
    
    if floorVal == ceilVal
        ceilVal = ceilVal + 1;
    end % END IF
    
    if ~(dataParams.classes(ii) == 11)
    
        trainIdx = trialIdx(randIdx(1:floorVal));
        testIdx = trialIdx(randIdx(ceilVal:end));
        
        trainTrial{ii} = logicalMat;
        testTrial{ii}  = logicalMat;
        
        trainTrial{ii}(trainIdx) = 1;
        testTrial{ii}(testIdx) = 1;
    
    end
end % END FOR uniqueTrials

idx11 = 0;
% Generate a list of approved 'silence' locations
if sum(dataParams.classes == 11)
    idx11 = find(dataParams.classes == 11);
    idxOther = find(dataParams.classes ~= 11,1);
    trainTrial{idx11} = trainTrial{idxOther};
    testTrial{idx11} = testTrial{idxOther};
    for ii = 1:size(dataParams.classes,2)
        if ~(ii==idx11)
            trainTrial{idx11} = trainTrial{idx11} | trainTrial{ii};
            testTrial{idx11}  = testTrial{idx11} | testTrial{ii};
        end
    end
end


%% Select and mold data

for jj = 1:size(dataParams.classes,2)
    if ~(jj == idx11)
        tmpDataP = squeeze(p(dataParams.timeBounds(1):dataParams.timeBounds(2),:,:,trainTrial{jj}));
        tmpDataP = permute(tmpDataP, [2,1,3]);
        tmpDataP = reshape(tmpDataP, size(tmpDataP,1), (diff(dataParams.timeBounds)+1)*size(tmpDataP,3));
        
        trainDataP{jj}  = tmpDataP;
        
        tmpDataP = squeeze(p(dataParams.timeBounds(1):dataParams.timeBounds(2),:,:,testTrial{jj}));
        tmpDataP = permute(tmpDataP, [2,1,3]);
        tmpDataP = reshape(tmpDataP, size(tmpDataP,1), (diff(dataParams.timeBounds)+1)*size(tmpDataP,3));
        
        testDataP{jj}  = tmpDataP;
    else
        tmpDataP = squeeze(p(1:8,:,:,trainTrial{jj}));
        tmpDataP = permute(tmpDataP, [2,1,3]);
        tmpDataP = reshape(tmpDataP, size(tmpDataP,1), 8*size(tmpDataP,3));
        
        trainZeroP{jj} = tmpDataP;
        
        tmpDataP = squeeze(p(1:8,:,:,testTrial{jj}));
        tmpDataP = permute(tmpDataP, [2,1,3]);
        tmpDataP = reshape(tmpDataP, size(tmpDataP,1), 8*size(tmpDataP,3));
        
        testZeroP{jj} = tmpDataP;
    end
end % END FOR uniqueTrials

sizeTrainClass = 0;

trainData  = [];
trainLabels   = [];

testData  = [];
testLabels   = [];

trainSilence = [];
testSilence = [];

for jj = 1:size(dataParams.classes,2)
    if ~(jj==idx11)
        trainData  = [trainData,  trainDataP{jj} ];
        trainLabels   = [trainLabels;  ones(size(trainDataP{jj},2),1) * dataParams.classes(jj)];
        testData  = [testData,  testDataP{jj} ];
        testLabels   = [testLabels;  ones(size(testDataP{jj},2),1) * dataParams.classes(jj) ];
    else
        trainSilence = [trainSilence, trainZeroP{jj}];
        testSilence = [testSilence, testZeroP{jj}];
    end
end % END FOR uniqueTrials

if idx11 >0
    % Assign Silence
    sizeTrainSilence = size(trainSilence,2);
    sizeTestSilence  = size(testSilence,2);
    
    sizeTrainClass  = size(trainDataP{1},2);
    sizeTestClass   = size(testDataP{1},2);
    
    availTrainSilence = randperm(sizeTrainSilence);
    availTestSilence  = randperm(sizeTestSilence);
    
    trainData = [trainData, trainSilence(:,availTrainSilence(1:sizeTrainClass))];
    testData  = [testData, testSilence(:,availTestSilence(1:sizeTestClass))];
    
    trainLabels = [trainLabels;  ones(sizeTrainClass,1) * dataParams.classes(idx11)];
    testLabels  = [testLabels;  ones(sizeTestClass,1) * dataParams.classes(idx11)];
end

% Randomly sort the rows and columns
testData = testData';
trainData = trainData';

rowSortTrain = [1:size(trainData,1)];
rowSortTest = [1:size(testData,1)];
colSort = [1:size(trainData,2)];

rowSortTrain = rowSortTrain(randperm(length(rowSortTrain)));
rowSortTest = rowSortTest(randperm(length(rowSortTest)));
colSort = colSort(randperm(length(colSort)));

trainLabels = trainLabels(rowSortTrain);
trainData = trainData(rowSortTrain,:);
trainData = trainData(:,colSort);

testLabels = testLabels(rowSortTest);
testData = testData(rowSortTest,:);
testData = testData(:,colSort);

trainData  = trainData.^2;

testData  = testData.^2;

testData = testData';
trainData = trainData';

end % END FUNCTION

% EOF