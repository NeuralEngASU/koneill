%% Delta Data


% fileList{1} = {'D:\PLI\Speech\DeltaSpeech_Day1_DS5000.mat'};
% fileList{2} = {'D:\PLI\Speech\DeltaSpeech_Day1_DS5000_BandPass250_2000.mat'};
% fileList{3} = {'D:\PLI\Speech\DeltaSpeech_Day1_DS5000_LowPass250Hz.mat'};
% fileList{1} = {'D:\PLI\Speech\DeltaSpeech_Day1_DS5000_LowPass250Hz_Notch60Hz.mat'};
% fileList{2} = {'D:\PLI\Speech\DeltaSpeech_Day1.mat'};
% fileList{3} = {'D:\PLI\Speech\DeltaSpeech_Day1_CombFiltered.mat'};

outputPath = 'D:\PLI\Speech\PLIOutput';

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
        [filePathOut] = GenPLISpeech(fileList{ii}{1}, outputPath, params);
    end % END FOR
end % END FOR

