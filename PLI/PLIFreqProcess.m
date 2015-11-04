%% Load raw data
szDataFile = '2014PP04Sz4.mat';
nonSzDataFile = '2014PP04NonSz4.mat';



% Load Raw Data
load(['E:\data\human CNS\EMD\Sz\clips\', szDataFile])
szData = detrend(data');

load(['E:\data\human CNS\EMD\NonSz\clips\', nonSzDataFile])
nonSzData = detrend(data');

clear data;

%% Low-Pass Filter to 250 Hz 

% Low pass 250 Hz Filter
Fs = 500;  % [Hz] Sampling Frequency

Fpass = 245;             % Passband Frequency
Fstop = 250;             % Stopband Frequency
Dpass = 0.057501127785;  % Passband Ripple
Dstop = 0.0001;          % Stopband Attenuation
dens  = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);
HdCoeffs = coeffs(Hd);

szDataFilt = zeros(size(szData));
nonSzDataFilt = zeros(size(nonSzData));

for ii = 1:size(szData,2)
    
    szDataFilt(:,ii)    = filtfilt(HdCoeffs.Numerator, 1, szData(:,ii));
    nonSzDataFilt(:,ii) = filtfilt(HdCoeffs.Numerator, 1, nonSzData(:,ii));
    disp(ii)
    
end % END FOR Low Pass Filter

%% Notch Filter 60Hz

Fs = 30000;
Fpass1 = 54;              % First Passband Frequency
Fstop1 = 59;              % First Stopband Frequency
Fstop2 = 61;              % Second Stopband Frequency
Fpass2 = 66;              % Second Passband Frequency
Dpass1 = 0.028774368332;  % First Passband Ripple
Dstop  = 0.001;           % Stopband Attenuation
Dpass2 = 0.057501127785;  % Second Passband Ripple
dens   = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass1 Fstop1 Fstop2 Fpass2]/(Fs/2), [1 0 ...
                          1], [Dpass1 Dstop Dpass2]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);

HdCoeffs = coeffs(Hd);


for ii = 1:size(szData,2)
    
    szDataFilt(:,ii)    = filtfilt(HdCoeffs.Numerator, 1, szDataFilt(:,ii));
    nonSzDataFilt(:,ii) = filtfilt(HdCoeffs.Numerator, 1, nonSzDataFilt(:,ii));
    disp(ii)
    
end % END FOR Low Pass Filter

%% Band Pass Filters

freqBandPass = [5,45; 55,95; 105,145; 155,195; 205,245];
freqBandStop = [freqBandPass(:,1)-5, freqBandPass(:,2)+5];

for jj = 2:5%:size(freqBandPass,1)
Fs = 500;  % Sampling Frequency

Fstop1 = freqBandStop(jj,1); % First Stopband Frequency
Fpass1 = freqBandPass(jj,1); % First Passband Frequency
Fpass2 = freqBandPass(jj,2); % Second Passband Frequency
Fstop2 = freqBandStop(jj,2); % Second Stopband Frequency
Dstop1 = 0.001;           % First Stopband Attenuation
Dpass  = 0.057501127785;  % Passband Ripple
Dstop2 = 0.0001;          % Second Stopband Attenuation
dens   = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 ...
                          0], [Dstop1 Dpass Dstop2]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);

HdCoeffs = coeffs(Hd);

eval(['szDataFilt', num2str(freqBandStop(jj,1)), '_', num2str(freqBandStop(jj,2)), '=zeros(size(szData));'])
eval(['nonSzDataFilt', num2str(freqBandStop(jj,1)), '_', num2str(freqBandStop(jj,2)), '=zeros(size(nonSzData));'])

for ii = 1:size(szData,2)
    
    
    eval(['szDataFilt', num2str(freqBandStop(jj,1)), '_', num2str(freqBandStop(jj,2)), '(:,ii)= filtfilt(HdCoeffs.Numerator, 1, szDataFilt(:,ii));'])
    eval(['nonSzDataFilt', num2str(freqBandStop(jj,1)), '_', num2str(freqBandStop(jj,2)), '(:,ii)= filtfilt(HdCoeffs.Numerator, 1, nonSzDataFilt(:,ii));'])

    disp(ii)
    
end % END FOR Low Pass Filter

szOutputFileName = ['D:\PLI\SeizureDetection\Processed\Sz\Clips\', szDataFile(1:end-4), '_',num2str(freqBandStop(jj,1)), '_', num2str(freqBandStop(jj,2)), '.mat'];
nonSzOtputFileName = ['D:\PLI\SeizureDetection\Processed\NonSz\Clips\', nonSzDataFile(1:end-4), '_',num2str(freqBandStop(jj,1)), '_', num2str(freqBandStop(jj,2)), '.mat'];

eval(['data = szDataFilt', num2str(freqBandStop(jj,1)), '_', num2str(freqBandStop(jj,2)), ';']);
save(szOutputFileName, 'data');
clear data

eval(['data = nonSzDataFilt', num2str(freqBandStop(jj,1)), '_', num2str(freqBandStop(jj,2)), ';']);
save(nonSzOtputFileName, ['nonSzDataFilt', num2str(freqBandStop(jj,1)), '_', num2str(freqBandStop(jj,2))]);
clear data

end % END FOR num freq windows

szOutputFileName = ['D:\PLI\SeizureDetection\Processed\Sz\Clips\', szDataFile(1:end-4), '_', 'Filt.mat'];
nonSzOtputFileName = ['D:\PLI\SeizureDetection\Processed\NonSz\Clips\', nonSzDataFile(1:end-4), '_', 'Filt.mat'];

data = szDataFilt;
save(szOutputFileName, 'data');
clear data

data = nonSzDataFilt;
save(nonSzOtputFileName, 'data');
clear data


%% Plot tests

figure
subplot(6,1,1)
plot(nonSzDataFilt(:,1))

subplot(6,1,2)
plot(nonSzDataFilt0_50(:,1))

subplot(6,1,3)
plot(nonSzDataFilt50_100(:,1))

subplot(6,1,4)
plot(nonSzDataFilt100_150(:,1))

subplot(6,1,5)
plot(nonSzDataFilt150_200(:,1))

subplot(6,1,6)
plot(nonSzDataFilt200_250(:,1))


%% PLI

fileList{1} = {'D:\PLI\SeizureDetection\Processed\Sz\Clips\2014PP04Sz4_100_150.mat'};
% fileList{2} = {'D:\PLI\SeizureDetection\Processed\Sz\Clips\2014PP04Sz4_50_100.mat'};
% fileList{3} = {'D:\PLI\SeizureDetection\Processed\Sz\Clips\2014PP04Sz4_150_200.mat'};
% fileList{4} = {'D:\PLI\SeizureDetection\Processed\Sz\Clips\2014PP04Sz4_200_250.mat'};
% fileList{5} = {'D:\PLI\SeizureDetection\Processed\Sz\Clips\2014PP04Sz4_Filt.mat'};
% fileList{6} = {'D:\PLI\SeizureDetection\Processed\NonSz\Clips\2014PP04NonSz4_0_50.mat'};
% fileList{7} = {'D:\PLI\SeizureDetection\Processed\NonSz\Clips\2014PP04NonSz4_50_100.mat'};
% fileList{8} = {'D:\PLI\SeizureDetection\Processed\NonSz\Clips\2014PP04NonSz4_150_200.mat'};
% fileList{9} = {'D:\PLI\SeizureDetection\Processed\NonSz\Clips\2014PP04NonSz4_200_250.mat'};
% fileList{10} = {'D:\PLI\SeizureDetection\Processed\NonSz\Clips\2014PP04NonSz4_Filt.mat'};


outputPath = 'D:\PLI\SeizureDetection\ProcessedPLI\';

winSize = 1;
Fs = 500;

for ii = 1:length(fileList)
    
    params.winSize = 1;
    params.Fs = 500;
    params.chanProcess = [];
    params.surrFlag = 0;
    params.surrNum = 0;
    params.rawPhiFlag = 0;
    params.biPolarFlag = 0;
    params.statsFlag = 0;
    params.globalFlag = 0;
    params.globalChan = [];

    [filePathOut] = GenPLIECoG(fileList{ii}{1}, outputPath, params);

end % END FOR
