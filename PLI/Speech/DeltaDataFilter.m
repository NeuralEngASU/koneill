%% Load Data

load('D:\PLI\Speech\DeltaSpeech_Day1.mat');
data = double(data);
data = detrend(data');

%% Design Filter

d = fdesign.comb('notch','L,BW,GBW,Nsh',10,5,-4,4,500);
Hd=design(d);
fvtool(Hd)


%% Apply Filter ??????

for ii = 1:32
    tic
    data(:,ii) = filtfilt(Hd.Numerator, Hd.Denominator, data3(:,ii));
    toc
end % END FOR

save('D:\PLI\Speech\DeltaSpeech_Day1_CombFiltered.mat', 'data', 'Header', 'lbl_words', 'offset_words', 'pts_words', '-v7.3');


%%

figure
subplot(2,1,1)
plot(data(23000000:24000000,1))
subplot(2,1,2)
plot(data2(23000000:24000000,1))

%% Downsample


load('D:\PLI\Speech\DeltaSpeech_Day1.mat');
data = double(data);
data = detrend(data');
data3 = data;

%% Low pass filter to 2000Hz

Fs = 30000;  % Sampling Frequency

Fpass = 2000;        % Passband Frequency
Fstop = 2400;        % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 80;          % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
Hd = design(h, 'butter', 'MatchExactly', match);

%% Filter data

for ii = 1:32
    tic
    data(:,ii) = filtfilt(Hd.sosMatrix, Hd.scaleValues, data3(:,ii));
    toc
    disp(ii)
end % END IF

%% Downsample to 5000Hz

Fs = 30000;
dsFs = 4800;

data = data(1:round(Fs/dsFs):end,:);

save('D:\PLI\Speech\DeltaSpeech_Day1_DS5000.mat', 'data', 'Header', 'lbl_words', 'offset_words', 'pts_words', '-v7.3');

data3 = data;
%% High Band Pass 250-2000

% All frequency values are in Hz.
Fs = 5000;  % Sampling Frequency

Fstop1 = 200;         % First Stopband Frequency
Fpass1 = 250;         % First Passband Frequency
Fpass2 = 2000;        % Second Passband Frequency
Fstop2 = 2050;        % Second Stopband Frequency
Astop1 = 60;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 80;          % Second Stopband Attenuation (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2, Fs);
Hd = design(h, 'butter', 'MatchExactly', match);


for ii = 1:32
    tic
    data(:,ii) = filtfilt(Hd.sosMatrix, Hd.scaleValues, data3(:,ii));
    toc
    disp(ii)
end % END IF
    
save('D:\PLI\Speech\DeltaSpeech_Day1_DS5000_BandPass250_2000.mat', 'data', 'Header', 'lbl_words', 'offset_words', 'pts_words', '-v7.3');


%% Low Pass

Fs = 5000;  % Sampling Frequency

Fpass = 250;         % Passband Frequency
Fstop = 300;         % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 80;          % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
Hd = design(h, 'butter', 'MatchExactly', match);

for ii = 1:32
    tic
    data(:,ii) = filtfilt(Hd.sosMatrix, Hd.scaleValues, data3(:,ii));
    toc
    disp(ii)
end % END IF
    
save('D:\PLI\Speech\DeltaSpeech_Day1_DS5000_LowPass250Hz.mat', 'data', 'Header', 'lbl_words', 'offset_words', 'pts_words', '-v7.3');



%% Notch Filter [60Hz]

data3 = data;

% All frequency values are in Hz.
Fs = 5000;  % Sampling Frequency

Fpass1 = 56;          % First Passband Frequency
Fstop1 = 59;          % First Stopband Frequency
Fstop2 = 61;          % Second Stopband Frequency
Fpass2 = 64;          % Second Passband Frequency
Apass1 = 0.5;         % First Passband Ripple (dB)
Astop  = 60;          % Stopband Attenuation (dB)
Apass2 = 1;           % Second Passband Ripple (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
                      Apass2, Fs);
Hd = design(h, 'butter', 'MatchExactly', match);

for ii = 1:32
    tic
    data(:,ii) = filtfilt(Hd.sosMatrix, Hd.scaleValues, data3(:,ii));
    toc
    disp(ii)
end % END IF
    
save('D:\PLI\Speech\DeltaSpeech_Day1_DS5000_LowPass250Hz_Notch60Hz.mat', 'data', 'Header', 'lbl_words', 'offset_words', 'pts_words', '-v7.3');



% EOF