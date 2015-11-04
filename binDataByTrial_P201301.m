function BinData = binDataByTrial_P201301(NEVIdx)

% Version date: 20130608
% Author: Tyler Davis

%% Loading a parsing trial data
FileStruct = loadFileStruct_P201301;
BinData.NEVFile = FileStruct.Decode(NEVIdx).NEV;
BinData.DiscardedOffline = FileStruct.Decode(NEVIdx).DiscardedOffline;
BinData.DiscardedOnline = FileStruct.Decode(NEVIdx).DiscardedOnline;

[NEVPath,NEVName] = fileparts(BinData.NEVFile);
BaseName = cell2mat(regexp(NEVName,'\d+-\d+-\d+','match'));
[TrainDS,RateDS] = calcDSData(BinData.NEVFile);
load(fullfile(NEVPath,[NEVName,'.mat']))
load(fullfile(NEVPath,[BaseName,'.ns5mat']),'-mat','Header')

TDataDS = TrainDS.Data/1000; %normalize to 1
RateDS = RateDS.Data/50;
TrialStruct = parseFingerPressPNS(NEV);

DOFList = {'ThumbFlex','IndexFlex','MiddleFlex','RingFlex','LittleFlex','IndexAbd','RingAbd','LittleAbd','ThumbExt','IndexExt','MiddleExt','RingExt','LittleExt'};

% Indexing offline trials. Restrict trials to training sets of 20 near beginning of
% recording.
Idx(1,:) = ~cellfun(@isempty,regexp({TrialStruct.MSMSSrc},'Training'));
Idx(2,:) = ~cellfun(@isempty,regexp({TrialStruct.TrainingSrc},'Cosine'));
Idx(3,:) = ~cellfun(@isempty,regexp({TrialStruct.Training},'Off'));
Idx(4,:) = ~cellfun(@isempty,{TrialStruct.TrialTypeName});
OfflineTrials = TrialStruct(all(Idx));
OfflineTrials(BinData.DiscardedOffline) = [];
if (all([OfflineTrials.CosHold]==[OfflineTrials.CosHold]) && all([OfflineTrials.CosRise]==[OfflineTrials.CosRise])) %make sure all trials are the same length
    %     TrialSamples = ceil((OfflineTrials(1).CosRise+OfflineTrials(1).CosHold)*15); %sampling frequency of 15 S/sec
    TrialSamples = 30; %fixed at 2sec to capture majority of movement and allow for baseline data to be collected without overlap
else
    disp('Trial timings vary!!!')
end
BinData.OfflineTrials = OfflineTrials;

% Indexing online trials
Idx(1,:) = ~cellfun(@isempty,regexp({TrialStruct.MSMSSrc},'Training'));
Idx(2,:) = ~cellfun(@isempty,regexp({TrialStruct.TrainingSrc},'Cosine'));
Idx(3,:) = ~cellfun(@isempty,regexp({TrialStruct.Training},'On'));
Idx(4,:) = ~cellfun(@isempty,regexp({TrialStruct.MSMSSrc},'Kalman'));
Idx(5,:) = ~cellfun(@isempty,regexp({TrialStruct.TrainingSrc},'Zero'));
Idx(6,:) = ~cellfun(@isempty,regexp({TrialStruct.Training},'Off'));
Idx(7,:) = ~cellfun(@isempty,{TrialStruct.TrialTypeName});
OnlineTrials = TrialStruct(all(Idx([1:3,7],:)) | all(Idx([4:6,7],:)));
OnlineTrials(BinData.DiscardedOnline) = [];
BinData.OnlineTrials = OnlineTrials;
    
% Find the number of trials for every movement type (dof) and locate the maximum
[DOF,DOFIdxs] = intersect(DOFList,{OfflineTrials.TrialTypeName});
[~,SortIdxs] = sort(DOFIdxs);
DOF = DOF(SortIdxs);
TrialCount = zeros(length(DOF),1);
for k=1:length(DOF)
    TrialCount(k) = sum(strcmp(DOF(k),{OfflineTrials.TrialTypeName}));
end
BinData.TrialCount = TrialCount;
BinData.DOF = DOF;
MaxTrials = max(TrialCount);

% Organize rates (features) and training variables (kinematics) into trial snippets for classification
Features = nan(length(DOF),96,MaxTrials,TrialSamples);
FeatureID = nan(length(DOF),96,MaxTrials,TrialSamples);
Baselines = nan(length(DOF),96,MaxTrials,TrialSamples);
BaselineASD = nan(length(DOF),96,MaxTrials,TrialSamples);
FKinematics = nan(length(DOF),MaxTrials,TrialSamples);
BKinematics = nan(length(DOF),MaxTrials,TrialSamples);
for k=1:length(DOF)
    DOFTrials = OfflineTrials(strcmp({OfflineTrials.TrialTypeName},DOF(k)));
    OnTS = ceil(double([DOFTrials.OnTS])/2000)'; %15 S/sec
  
    Offset = 7; %500ms offset from target appearance to capture start of movement
    MTS = repmat(0:(TrialSamples-1),length(OnTS),1) + repmat(OnTS,1,TrialSamples) + Offset; %timestamp matrix for movement period    
    BTS = repmat((1-TrialSamples):0,length(OnTS),1) + repmat(OnTS,1,TrialSamples); %timestamp matrix for baseline period
    
    Features(k,:,1:size(MTS,1),:) = reshape(RateDS(:,MTS),[96,size(MTS)]); %data snippets of length "TrialSamples" with offset of 500ms from target appearance (dof x electrode x trial x rate)
    FeatureID(k,:,1:size(MTS,1),:) = zeros(96,size(MTS,1),TrialSamples)+k; %index into DOF cell array
    Baselines(k,:,1:size(BTS,1),:) = reshape(RateDS(:,BTS),[96,size(BTS)]); %data snippets of length "TrialSamples" before movement (dof x electrode x trial x rate)
    BaselineASD(k,:,1:size(BTS,1),:) = repmat([DOFTrials.ASD],[96,1,TrialSamples]); %acquire start durations for all trials
    if k>=9
        FKinematics(k,1:size(MTS,1),:) = reshape(TDataDS(k-8,MTS),size(MTS)); %kinematic snippets that correspond with features (dof x trial x rate)
        BKinematics(k,1:size(BTS,1),:) = reshape(TDataDS(k-8,BTS),size(BTS)); %kinematic snippets that corrspond with baselines (dof x trial x rate)
    else       
        FKinematics(k,1:size(MTS,1),:) = reshape(TDataDS(k,MTS),size(MTS));
        BKinematics(k,1:size(BTS,1),:) = reshape(TDataDS(k,BTS),size(BTS));
    end
end
BinData.Features = Features;
BinData.FeatureID = FeatureID;
BinData.Baselines = Baselines;
BinData.BaselineASD = BaselineASD;
BinData.FKinematics = FKinematics;
BinData.BKinematics = BKinematics;
if isfield(FileStruct.Decode,'DrivenElects')
    if ~isempty(FileStruct.Decode(NEVIdx).DrivenElects)
        BinData.DrivenElectrodes = FileStruct.Decode(NEVIdx).DrivenElects;
    else
        BinData.DrivenElectrodes = findDrivenElects(Features,Baselines);
    end
else
    BinData.DrivenElectrodes = findDrivenElects(Features,Baselines);
end

% Kalman offline parsing
KTrainIdxs = cell(1,length(DOF));
KTestIdxs = cell(1,length(DOF));
KTrainFeatures = cell(1,length(DOF));
KTestFeatures = cell(1,length(DOF));
KTrainKinematics = cell(1,length(DOF));
KTestKinematics = cell(1,length(DOF));
for k=1:length(DOF)
    DOFTrials = OfflineTrials(strcmp({OfflineTrials.TrialTypeName},DOF(k)));
    OnTS = ceil(double([DOFTrials.OnTS])/2000)'; %15 S/sec
    
    KTrainIdxs(k) = {OnTS(1):OnTS(10)+75};
    KTestIdxs(k) = {OnTS(11):OnTS(end)+75};
    KTrainKinematics(k) = {zeros(length(DOF),length(KTrainIdxs{k}))};
    KTestKinematics(k) = {zeros(length(DOF),length(KTestIdxs{k}))};
    
    KTrainFeatures(k) = {RateDS(:,KTrainIdxs{k})};
    KTestFeatures(k) = {RateDS(:,KTestIdxs{k})};
    if k>=9
        KTrainKinematics{k}(k,:) = TDataDS(k-8,KTrainIdxs{k});
        KTestKinematics{k}(k,:) = TDataDS(k-8,KTestIdxs{k});
    else
        KTrainKinematics{k}(k,:) = TDataDS(k,KTrainIdxs{k});
        KTestKinematics{k}(k,:) = TDataDS(k,KTestIdxs{k});
    end
end

% Adding a rest period as 14th DOF. This rest period is a concatenation of
% a 1sec window immediately preceding the OnTS marker for the 1st 6 trials
% of the 13 unique movements. This rest is only added to the testing set.
KTestFeatures = [KTestFeatures,{reshape(permute(Baselines(:,:,1:6,15:end),[2,4,1,3]),96,[])}];
KTestKinematics = [KTestKinematics,{zeros(length(DOF),size(KTestFeatures{end},2))}];

BinData.KTrainIdxs = KTrainIdxs;
BinData.KTestIdxs = KTestIdxs;
BinData.KTrainFeatures = KTrainFeatures;
BinData.KTestFeatures = KTestFeatures;
BinData.KTrainKinematics = KTrainKinematics;
BinData.KTestKinematics = KTestKinematics;

% Kalman online
if ~isempty(BinData.DiscardedOnline)
    OnTS = ceil(double([OnlineTrials.OnTS])/2000); %15 S/sec
    OnlineIdxs = (OnTS(1):OnTS(end)+75);
    BinData.KOnlineKinematics = TDataDS(:,OnlineIdxs);
end



%% Plotting trial-binned data for verification
% clf;
% trial = 4; mvnt = 3; elect = 44;
% subplot(1,2,1);
% hold on;
% plot(squeeze(Baselines(mvnt,e2c(elect,'pns'),trial,:)));
% plot(squeeze(BKinematics(mvnt,trial,:)));
% hold off;
% axis([1,TrialSamples,-1,1]);
% title([DOF{mvnt},' (',num2str(squeeze(BaselineASD(mvnt,e2c(elect,'pns'),trial,1))),')'])
% 
% subplot(1,2,2);
% hold on;
% plot(squeeze(Features(mvnt,e2c(elect,'pns'),trial,:)));
% plot(squeeze(FKinematics(mvnt,trial,:)));
% hold off;
% axis([1,TrialSamples,-1,1])
% title(DOF{mvnt})





