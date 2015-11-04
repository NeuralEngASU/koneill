function corrOut = PNSCorr( filePath, trialStruct, sampleRate )

corrOut = 1;

% Open relevent files
nsx2mat([filePath,'.ns5']);
NEV = openNEV([filePath,'.nev']);
load([filePath,'.ns5mat'],'-mat','Header')

% Down samples data.
sampleSkip = 1;

% Detect if sample rate is too high
if sampleRate > Header.Fs
    error('Sample rate is too high, must be less than Fs of data')
elseif sampleRate < 60
    error('Sample rate is too low, must be greater than 60')
else
    sampleSkip = ceil(Header.Fs/sampleRate);
    disp(sprintf('Analyzing data at %0.2f S/sec', Header.Fs/sampleSkip))
end % END IF

% Initialize DOFs
idxDOF = 137:144;
DOF = {'Thumb','Index','Middle','Ring','Little','IndexIntrinsic','RingIntrinsic','LittleIntrinsic'};

%
if exist([filePath,'_TDataDS.mat'],'file')
    load([filePath,'_TDataDS.mat'])
else
    Hd = design(fdesign.lowpass('N,F3dB',4,15,30000),'butter');
    TDataDS = zeros(8,ceil(Header.ChannelLengthSamples/sampleSkip)); % Total length when downsampled
    for k=1:length(idxDOF)
        TData = readNSxMat([filePath,'.ns5mat'],idxDOF(k));
        TData = filter(Hd,TData.data);
        TDataDS(k,:) = TData(1:sampleSkip:end); % Downsample to inputed value
    end
    clear('TData')        
    save([filePath,'_TDataDS.mat'],'TDataDS')
end
TDataDS = TDataDS/1000; %normalize to 1


%%%%% I do not know why the RateB veriable is down sampled to 300 S/sec
%%%%% for the time being, use 60 S/sec as target fs
if exist([filePath,'_RateDS.mat'],'file')
    load([filePath,'_RateDS.mat'])
elseif sampleRate == 60
    boxwin = 0.3; % boxcar window in sec
    TSCell = cell(96,1);
    WFCell = cell(96,1);
    RateDS = zeros(96,ceil(Header.ChannelLengthSamples/sampleSkip)); % Total length when downsampled
    Hd = design(fdesign.lowpass('N,F3dB',4,15,300),'butter');
    for k=1:96
        clc, disp(k)
        TS = double(NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Electrode==k)); 
        TSCell(k) = {ceil(TS./sampleSkip)}; %timestamps at downsample input rate
        WFCell(k) = {double(NEV.Data.Spikes.Waveform(:,NEV.Data.Spikes.Electrode==k))};
        RateB = zeros(1,ceil(Header.ChannelLengthSamples/100)); %total length for 300 S/sec
        if ~isempty(TS)
            RateB(ceil(TS./100)) = 1; %timestamps at 300 S/sec
        end
        %     RateB = conv2(RateB,ones(1,boxwin*60)./boxwin,'same');
        RateB = conv2(RateB,gausswin(boxwin*300)'./(boxwin/2),'same'); %std = N/5 (i.e. std = 90/5 = 18samples = 60ms)
        RateB = filter(Hd,RateB); %filter at 15Hz
        RateDS(k,:) = RateB(1:5:end); %downsample from 300 S/sec to 60 S/sec
    end
    clear('RateB')        
    save([filePath,'_RateDS.mat'],'TSCell','WFCell','RateDS')
elseif sampleRate == 10000
    boxwin = 0.3; % boxcar window in sec
    TSCell = cell(96,1);
    WFCell = cell(96,1);
    RateDS = zeros(96,ceil(Header.ChannelLengthSamples/sampleSkip)); % Total length when downsampled
    Hd = design(fdesign.lowpass('N,F3dB',4,15,300),'butter');
    for k=1:96
        clc, disp(k)
        TS = double(NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Electrode==k)); 
        TSCell(k) = {ceil(TS./sampleSkip)}; %timestamps at downsample input rate
        WFCell(k) = {double(NEV.Data.Spikes.Waveform(:,NEV.Data.Spikes.Electrode==k))};
        RateB = zeros(1,ceil(Header.ChannelLengthSamples/sampleSkip)); %total length for 300 S/sec
        if ~isempty(TS)
            RateB(ceil(TS./100)) = 1; %timestamps at 300 S/sec
        end
        %     RateB = conv2(RateB,ones(1,boxwin*60)./boxwin,'same');
        RateB = conv2(RateB,gausswin(boxwin*300)'./(boxwin/2),'same'); %std = N/5 (i.e. std = 90/5 = 18samples = 60ms)
        RateB = filter(Hd,RateB); %filter at 15Hz
        RateDS(k,:) = RateB(1:1:end); %downsample from 300 S/sec to 60 S/sec
    end
    clear('RateB')        
    save([filePath,'_RateDS.mat'],'TSCell','WFCell','RateDS')
end
TSCell = TSCell;
WFCell = WFCell; %unit waveforms
RateDS = RateDS/50; %normalize so 50Hz is 1

fieldList = fieldnames(trialStruct);
fieldList2DOFMap = [1, 8, 7, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1];

corrOut = zeros(96, size(fieldList,1));

for k = 1:size(fieldList, 1)
%     if k==6
%         disp('Hello World')
%     end
    
    begTime = eval(['trialStruct.',fieldList{k},'(1,  4);']);
    endTime = eval(['trialStruct.',fieldList{k},'(20, 4);']);

    trialwin = (endTime - begTime)/sampleSkip; % Down sample to input rate
    plotIdxs = 1:trialwin + begTime/sampleSkip;
    
%     [~,currdof] = max(abs(mean(TDataDS(:,plotidxs),2)));
    currDOF = fieldList2DOFMap(k);
    
    mTDataDS = TDataDS(currDOF,plotIdxs) - repmat(mean(TDataDS(currDOF,plotIdxs),2),1,length(plotIdxs));
    sTDataDS = std(TDataDS(currDOF,plotIdxs),0,2);
    
    mRateDS = RateDS(:,plotIdxs) - repmat(mean(RateDS(:,plotIdxs),2),1,length(plotIdxs));
    sRateDS = std(RateDS(:,plotIdxs),0,2);
    
    C = abs((mRateDS*mTDataDS')/(length(plotIdxs)-1)); %covariance
    
    CC = C./(sRateDS*sTDataDS'); %pearson's correlation coefficient
    CC(isnan(CC))=0; CC(isinf(CC))=0;
    
    corrOut(:, k) = CC;

end % END FOR
end % END FUNCTION

% EOF