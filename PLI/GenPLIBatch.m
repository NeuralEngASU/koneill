
% fileList{1} = {'E:\data\human CNS\EMD\Sz\clips\2014PP02Sz1.mat'};
% fileList{2} = {'E:\data\human CNS\EMD\Sz\clips\2014PP02Sz2.mat'};
% fileList{3} = {'E:\data\human CNS\EMD\Sz\clips\2014PP02Sz3.mat'};
% fileList{4} = {'E:\data\human CNS\EMD\Sz\clips\2014PP02Sz4.mat'};
% fileList{5} = {'E:\data\human CNS\EMD\Sz\clips\2014PP02Sz5.mat'};
% fileList{6} = {'E:\data\human CNS\EMD\Sz\clips\2014PP02Sz6.mat'};

% fileList{1} = {'E:\data\human CNS\EMD\NonSz\ProcData\DN\2014PP01NonSz1_DN.mat'};
% fileList{2} = {'E:\data\human CNS\EMD\NonSz\ProcData\DN\2014PP01NonSz2_DN.mat'};
% fileList{3} = {'E:\data\human CNS\EMD\NonSz\ProcData\DN\2014PP01NonSz3_DN.mat'};
% fileList{1} = {'E:\data\human CNS\EMD\Sz\ProcData\DN\2014PP01Sz7_DN.mat'};
% fileList{1} = {'E:\data\human CNS\PLI_long_data\LongPLIclip1.mat'};
% fileList{1} = {'E:\data\human CNS\PLI_long_data\LongPLIclip3.mat'};
% fileList{3} = {'E:\data\human CNS\EMD\Sz\ProcData\DN\2014PP01Sz3_DN.mat'};

outputPath = 'E:\data\PLI\EMDData\LongFormData';

winSize = [1];
Fs = [500,500, 500, 500, 500, 500, 500];

for ii = 1:length(fileList)
    for jj = 1:length(winSize)
        [filePathOut] = GenPLI(fileList{ii}{1}, outputPath, 'winSize', winSize(jj), 'GlobalFlag', 1, 'rawPhiFlag',0, 'Fs', Fs(ii), 'STATSFLAG', 0);
    end % END FOR
end % END FOR

%% Delta Data Freq

%% Delta Data

fileList{1} = {'E:\data\PLI\delta\Verbal\20080730-151259\20080730-151259-002.ns5'};

outputPath = 'E:\data\PLI\delta\PLIOutput';

params.winSize = 0.0125;
params.Fs = 30000;
params.chanProcess = [1:32];
params.surrFlag = 0;
params.surrNum = 100;
params.rawPhiFlag = 0;
params.biPolarFlag = 0;
params.statsFlag = 0;
params.globalFlag = 0;
params.globalChan = [1:32];
params.word = 'yes';

for kk = 1:10

    wordList = {'yes'    'no'    'hot'    'cold'    'hungry'    'thirsty'    'hello'    'goodbye'    'more'    'less'};
    
    params.word = wordList{kk};
for ii = 1:length(fileList)
    for jj = 1:length(params.winSize)
        [filePathOut] = GenPLIVerbal(fileList{ii}{1}, outputPath, params);
    end % END FOR
end % END FOR
end % END FOR WORD
%%
params.winSize = 1;
params.Fs = 500;
params.chanProcess = [1:32];
params.surrFlag = 0;
params.surrNum = 0;
params.rawPhiFlag = 0;
params.biPolarFlag = 0;
params.statsFlag = 0;
params.globalFlag = 0;
params.globalChan = [1:32];

for ii = 1:length(fileList)
    for jj = 1:length(params.winSize)
        [filePathOut] = GenPLIChan(fileList{ii}{1}, outputPath, params);
    end % END FOR
end % END FOR

params.winSize = 0.25;
params.Fs = 30000;
params.chanProcess = [1:32];
params.surrFlag = 0;
params.surrNum = 100;
params.rawPhiFlag = 0;
params.biPolarFlag = 0;
params.statsFlag = 0;
params.globalFlag = 0;
params.globalChan = [1:32];

for ii = 1:length(fileList)
    for jj = 1:length(params.winSize)
        [filePathOut] = GenPLIChan(fileList{ii}{1}, outputPath, params);
    end % END FOR
end % END FOR

%% Delta Bands + Original

fileList{1} = {'D:\PLI\Delta_ProcessedTrialData.mat'};
% fileList{1} = {'E:\data\PLI\delta\PLIOutput\SeperateBands\Delta_PLI_Alpha.mat'};
% fileList{1} = {'E:\data\PLI\delta\PLIOutput\SeperateBands\Delta_PLI_Beta.mat'};
% fileList{1} = {'E:\data\PLI\delta\PLIOutput\SeperateBands\Delta_PLI_Chi.mat'};
% fileList{2} = {'E:\data\PLI\delta\PLIOutput\SeperateBands\Delta_PLI_Delta.mat'};
% fileList{3} = {'E:\data\PLI\delta\PLIOutput\SeperateBands\Delta_PLI_Gamma.mat'};
% fileList{4} = {'E:\data\PLI\delta\PLIOutput\SeperateBands\Delta_PLI_Theta.mat'};

outputPath = 'D:\PLI';

params.winSize = 0.1;
params.Fs = 5000;
params.chanProcess = [1:32];
params.surrFlag = 0;
params.surrNum = 0;
params.rawPhiFlag = 0;
params.biPolarFlag = 0;
params.statsFlag = 0;
params.globalFlag = 0;
params.globalChan = [1:32];

for ii = 1:length(fileList)
    for jj = 1:length(params.winSize)
        [filePathOut] = GenPLIVerbal(fileList{ii}{1}, outputPath, params);
    end % END FOR
end % END FOR

%% Batch Process Seizure PLI

szFilePath = 'E:\data\human CNS\EMD\Sz\clips';
szOutputPath = 'D:\PLI\SeizureDetection\Sz\HilbertFirst';
nonSzFilePath = 'E:\data\human CNS\EMD\NonSz\clips';
nonSzOutputPath = 'D:\PLI\SeizureDetection\NonSz\HilbertFirst';

szFileName = dir([szFilePath, '\*.mat']);
nonSzFileName = dir([nonSzFilePath, '\*.mat']);
numFiles = size(szFileName,1);

for ii =15:numFiles % SKIP NUMBER 15 IN CLIPS
    
    load(fullfile(szFilePath, szFileName(ii).name));
    numChans = size(data,1);
    
    params.winSize = 1;
    params.Fs = 500;
    params.chanProcess = [1:numChans];
    params.surrFlag = 0;
    params.surrNum = 0;
    params.rawPhiFlag = 0;
    params.biPolarFlag = 0;
    params.statsFlag = 0;
    params.globalFlag = 0;
    params.globalChan = [1:numChans];
    
    [~] = GenPLIECoG(fullfile(   szFilePath,    szFileName(ii).name),    szOutputPath, params);
    
    load(fullfile(nonSzFilePath, nonSzFileName(ii).name));
    numChans = size(data,1);
    
    params.winSize = 1;
    params.Fs = 500;
    params.chanProcess = [1:numChans];
    params.surrFlag = 0;
    params.surrNum = 0;
    params.rawPhiFlag = 0;
    params.biPolarFlag = 0;
    params.statsFlag = 0;
    params.globalFlag = 0;
    params.globalChan = [1:numChans];
    
    [~] = GenPLIECoG(fullfile(nonSzFilePath, nonSzFileName(ii).name), nonSzOutputPath, params);
    
    
end % END FOR


% EOF