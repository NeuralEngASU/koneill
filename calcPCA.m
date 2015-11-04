function PCA = calcPCA(NEVIdx,ClassIdxs)

% Version date: 20130608
% Author: Tyler Davis
%
% A set of indices (ClassIdxs, ranging from 1-14) are chosen for
% classification, which correspond with the 13 unique movements and rest
% given in ClassList. Training sets always include the first 10 trials for
% each movement type specified in ClassIdxs. Testing sets include the
% remaining number of trials for the same set of movement types, which is
% usually 10 plus or minus a few trials. ClassIdxs may or may not include
% rest (index=14), which is always placed at the end of the set. If rest is
% included, the baseline (rest) parsing code is executed. Rest periods are a 2sec
% window of time preceding the presentation of the target that indicates
% movement (OnTS). If more than one movement type is included in the set
% (i.e. thumbflex, indexflex, etc.), then there are more than 10 baseline
% trials to use for both the training and testing sets. Under this
% condition, random indices are generated to use for selecting 10 baseline
% trials for training and 10 baseline trials for testing. Movement periods
% are a 2sec window of time 500ms after the presentation of the target
% (OnTS). Movement periods correspond to the 13 different movement types
% given in ClassList (minus "Rest"). For the different movements, similar
% to the baselines, the first 10 trials are used for training and the last
% 10 trials for testing. The testing trials may vary from 10 by a few
% trials since exactly 20 trials total were not performed during data
% collection. To create the training and testing matrices, data is parsed
% by electrode, trial, and movement type and shaped into a matrix where the
% rows are a vertical concatenation of 10 trials for each of the different
% classes and the columns are a horizontal concatenation of snippets (2sec) of
% firing rates for each of the specified electrodes. Electrodes are chosen
% that show a significant increase in rate when compared with baseline. PCA
% is then performed on the training matrix to create a transformation
% matrix that is used to project both the training and testing data into PC
% space. This rotates these data sets so that the most variation between
% trials is represented by the first few dimensions (columns) of the
% transformed matrices. The first several columns that represent at most 90
% percent of the variation of the training set are used to classify the
% testing set. Classification is performed using a linear discriminate
% method that fits a multivariate normal density to each group, with a
% pooled estimate of covariance (seeks to project data to a line in
% n-dimensional space that maximizes the separation of group means, but
% minimizes group scatter or variance).

try
    %% Binning data by trial for PCA analysis
    BinData = binDataByTrial(NEVIdx); %NEVIdx corresponds to a specific dataset that is hardcoded into "binDataByTrial"   
    Baselines = BinData.Baselines;
    Features = BinData.Features;
    FeatureID = BinData.FeatureID;
    Elects = BinData.DrivenElectrodes;
    % Elects = Elects(1:20);
    % Elects = [6,13,36,41,44,45,79,88,89];
%     Elects = [2:10,12:80,82:90,92:100];
    ClassList = [BinData.DOF,{'Rest'}];   
    
    %% Selecting electrodes and movement types and performing PCA-based classification
    RStruct = struct('TrainMat',[],'TrainMatID',[],'TestMat',[],'TestMatID',[],'PCTrainMat',[],'PCTestMat',[],'TestResults',[],'PCorrect',[],'ClassIdxs',[],'Warning',[]);
    switch nargin
        case 1
            matlabpool 8
            TSize = zeros(1,length(ClassList)-1); for k=1:length(ClassList)-1; TSize(k) = nchoosek(length(ClassList),k+1); end;
            PTotals = nan(max(TSize),length(ClassList)-1);
            RStruct = repmat(RStruct,max(TSize),length(ClassList)-1);
            for m=1:length(ClassList)-1
                clc; disp(m)
                DOFIdxsMat = nchoosek(1:length(ClassList),m+1);
                for k=1:size(DOFIdxsMat,1)
                    ClassIdxs = DOFIdxsMat(k,:);
                    Results = runClassification(ClassIdxs,Elects,Baselines,Features,FeatureID);
                    PTotals(k,m) = Results.PCorrect;
                    if length(ClassIdxs)==length(ClassList)
                        RStruct(k,m) = Results;
                    end
                end
            end
            matlabpool close
        case 2
            Results = runClassification(ClassIdxs,Elects,Baselines,Features,FeatureID);
            PTotals = Results.PCorrect;
            RStruct = Results;
    end
    
    PCA.Results = RStruct(1,end);
    PCA.Electrodes = Elects;
    PCA.PTotals = PTotals;
    PCA.Trials = BinData.OfflineTrials;
    PCA.TrialCount = BinData.TrialCount;
    PCA.NEVFile = BinData.NEVFile;
    PCA.ClassList = ClassList;
    PCA.Baselines = Baselines;
    PCA.Features = Features;
    PCA.FeatureID = FeatureID;
    PCA.FKinematics = BinData.FKinematics;
    PCA.BKinematics = BinData.BKinematics;
    [NEVPath,NEVName] = fileparts(BinData.NEVFile);
    if nargin==1
        save(fullfile(NEVPath,['PCA_',NEVName,'_',datestr(clock,'yyyymmdd-HHMMSS'),'.mat']),'PCA')
    end
catch ME
    assignin('base','ME',ME)    
end


%%%%%%%%%%%%%%%%%%%%%% Subfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Results = runClassification(ClassIdxs,Elects,Baselines,Features,FeatureID)

RestIdx = size(Features,1)+1; %index value for rest periods (14 with finger extensions; 9 without extensions)
FeatureIdxs = ClassIdxs(ClassIdxs~=RestIdx);
%%%%%%%%%%%%%%%%%%%%%%% Randomly choosing basline data %%%%%%%%%%%%%%%%%%%
RTrain = []; RTrainID = []; RTest = []; RTestID = [];
if length(ClassIdxs)~=length(FeatureIdxs)
    RTrain = reshape(permute(Baselines(FeatureIdxs,e2c(Elects,'pns'),1:10,:),[3,1,4,2]),10*length(FeatureIdxs),[]);
    RTest = reshape(permute(Baselines(FeatureIdxs,e2c(Elects,'pns'),11:20,:),[3,1,4,2]),10*length(FeatureIdxs),[]);
    if length(FeatureIdxs)==1
        RTest = RTest(~all(isnan(RTest),2),:);
    else
        RTrain = RTrain(~all(isnan(RTest),2),:); RTest = RTest(~all(isnan(RTest),2),:);
        RIdxs = randperm(size(RTrain,1));
        RTrain = RTrain(RIdxs(1:10),:); RTest = RTest(RIdxs(1:10),:);
    end
    RTrainID = zeros(size(RTrain))+RestIdx; RTestID = zeros(size(RTest))+RestIdx;
    ClassIdxs = [FeatureIdxs,RestIdx];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Selecting movement data %%%%%%%%%%%%%%%%%%%%%%%%%%
FTrain = reshape(permute(Features(FeatureIdxs,e2c(Elects,'pns'),1:10,:),[3,1,4,2]),10*length(FeatureIdxs),[]); %dof x electrode x trial x rate -> trial(changes fastest down columns of desired final matrix) x dof(changes 2nd fastest down columns) x rate(etc) x electrode(changes slowest down columns); Final matrix is: trial+dof(blocks) x rate+electrode(blocks)
FTrainID = reshape(permute(FeatureID(FeatureIdxs,e2c(Elects,'pns'),1:10,:),[3,1,4,2]),10*length(FeatureIdxs),[]);
FTest = reshape(permute(Features(FeatureIdxs,e2c(Elects,'pns'),11:20,:),[3,1,4,2]),10*length(FeatureIdxs),[]);
FTestID = reshape(permute(FeatureID(FeatureIdxs,e2c(Elects,'pns'),11:20,:),[3,1,4,2]),10*length(FeatureIdxs),[]);
FTestID = FTestID(~all(isnan(FTest),2),:); FTest = FTest(~all(isnan(FTest),2),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TrainMat = [FTrain;RTrain];
TrainMatID = [FTrainID;RTrainID];
TestMat = [FTest;RTest];
TestMatID = [FTestID;RTestID];

%%%%%%%%%%%%%%%%%%%%% Shuffle testing data set %%%%%%%%%%%%
% TestMat = reshape(TestMat',TrialSamples,[]);
% TestMat = TestMat(:,randperm(size(TestMat,2)));
% TestMat = reshape(TestMat,length(Elects)*TrialSamples,[])';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[V,D] = eig(cov(TrainMat));
V = fliplr(V); D = fliplr(D); %organizing vectors so most variance is in left column
Var90 = find(cumsum(sum(D))/sum(D(:))<=0.9,1,'last'); %finding first set of vectors with 90 percent of the variance
if isempty(Var90)
    Var90 = 1;
end

PCTrainMat = TrainMat*V(:,1:Var90); %projecting to pc space
PCTestMat = TestMat*V(:,1:Var90);

w = '';
if any(sum(PCTrainMat,2)==0) %checking for null trials that have a zero value
    w = 'Insufficient dataset for classification';
    disp(w);
end

%%%%%%%%%%%%%%%% Nearest Centroid Classification %%%%%%%%%%%%%%%%%
% a=repmat(squeeze(mean(reshape(PCTrainMat,10,length(ClassIdxs),Var90))),[1,1,size(PCTestMat,1)]);
% b=permute(repmat(PCTestMat,[1,1,length(ClassIdxs)]),[3,2,1]);
% [~,C]=min(sqrt(squeeze(sum((a-b).^2,2))));
% C_Results = ClassIdxs(C)'==TestMatID(:,1);
% C_Total(k,m) = sum(C_Results)/length(C);
% C_Total = sum(C_Results)/length(C);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% Linear Discriminant Classification %%%%%%%%%%%%%%%%%%
L = classify(PCTestMat,PCTrainMat,TrainMatID(:,1)); %classifying
L_Results = (L==TestMatID(:,1));
L_Total = sum(L_Results)/length(L); %calculating percentage correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Results.TrainMat = TrainMat;
Results.TrainMatID = TrainMatID;
Results.TestMat = TestMat;
Results.TestMatID = TestMatID;
Results.PCTrainMat = PCTrainMat;
Results.PCTestMat = PCTestMat;
Results.TestResults = L;
Results.PCorrect = L_Total;
Results.ClassIdxs = ClassIdxs;
Results.Warning = w;
