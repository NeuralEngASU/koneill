%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Epilepsy Frequency Process
%   Lab: Neural Engineering Lab
%   PI: Bradley Greger
%   Author: Kevin O'Neill
%   Date: 2015.11.12
%
%   Desc: Filter (and optional bandpass) data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [fileList] = EpilepsyFreqProcess(filePath, fileName, params, Header)
% Load Raw Data
disp('*******************************')
disp('Loading Data')

load(fullfile(filePath, fileName))
if size(data,1) > size(data, 2)
    data = detrend(data);
else
    data = detrend(data');
end % END IF

dataFilt = zeros(size(data));

Fs = Header.Fs;

fileList = [];

%% Comb Filter 60Hz, 120Hz, 180Hz, 240Hz,...
if params.combFilter % Run if the user has selected combFilter
    
    disp('*******************************')
    disp('Comb Filtering Data')
    
    if isempty(params.combFilterBand)
        combHz = [60:60:(Fs/2)];
    else
        combHz = params.combFilterBand;
    end % END IF select comb filter frequencies;
    
    totalTimeTic = tic;
    
    dataFilt = data;
    
    for notchIdx = 1:length(combHz)
        
        notchHz = combHz(notchIdx);
        
        Fpass1 = notchHz-6;   % First Passband Frequency
        Fstop1 = notchHz-1;   % First Stopband Frequency
        Fstop2 = notchHz+1;   % Second Stopband Frequency
        Fpass2 = notchHz+6;   % Second Passband Frequency
        Apass1 = 1;         % First Passband Ripple (dB)
        Astop  = 10;          % Stopband Attenuation (dB)
        Apass2 = 1;           % Second Passband Ripple (dB)
        match  = 'stopband';  % Band to match exactly
        
        % Construct an FDESIGN object and call its BUTTER method.
        h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
            Apass2, Fs);
        Hd = design(h, 'butter', 'MatchExactly', match);     
        
        for chan = 1:size(dataFilt,2)
            
            dataFilt(:,chan)    = filtfilt(Hd.sosMatrix, Hd.ScaleValues, dataFilt(:,chan));
            
            totalTime = toc(totalTimeTic);
            
            numChansComplete = chan + size(dataFilt,2)*(notchIdx-1);
            
            chanTime = totalTime/numChansComplete;
            timeLeft = chanTime * (size(dataFilt,2)-chan) + chanTime * size(dataFilt,2) * (length(combHz)-notchIdx);
            
            clc;
            disp('*******************************')
            disp('Comb Filtering Data')
            disp([num2str(notchHz),'Hz Notch Filter. Chan: ', num2str(chan)])
            disp(['Time per chan: ', num2str(chanTime), ' sec'])
            disp(['Time spent: ', num2str(totalTime), ' sec'])
            disp(['Time Left: ', num2str(timeLeft), ' sec'])
        end % END FOR Notch Channel
    end % END FOR Notch Filter
    
    Header.combParams.FiltType = 'butter';
    Header.combParams.Fstop = combHz;
    Header.combParams.Fpass1 = -6;
    Header.combParams.Fstop1 = -1;
    Header.combParams.Fstop2 = +1;
    Header.combParams.Fpass2 = +6;
    Header.combParams.Astop = Astop;
    Header.combParams.Apass1 = 1;
    Header.combParams.Apass2 = 1;    
    Header.combParams.match = 'stopband';
    
    combOutputFullpath = [params.targetDir, fileName(1:end-4), '_', 'Comb.mat'];

    Header.fileName = [fileName(1:end-4), '_', 'Comb.mat'];
    Header.filePath = params.targetDir;
    
    data = dataFilt;
    save(combOutputFullpath, 'data', 'Header');
    clear data
    
    fileList{end+1} = Header.fileName;
        
else
    dataFilt = data;
end % END IF params.combFilter

%% Band Pass Filters

if params.bandPassFilter % Run if the user has selected to bandpass the data
    
    freqBandStop = params.passBands;
    freqBandPass = [freqBandStop(:,1)+5, freqBandStop(:,2)-5];    
    
    totalTimeTic = tic;
    
    for band = 1:size(freqBandPass,1)
        
        Fstop1 = freqBandStop(band,1); % First Stopband Frequency
        
%         if Fstop1 == 0
%             Fstop1 = 0.01;
%         end % END IF
        
        Fpass1 = freqBandPass(band,1);   % First Passband Frequency
        Fpass2 = freqBandPass(band,2);   % Second Passband Frequency
        Fstop2 = freqBandStop(band,2);   % Second Stopband Frequency
        Astop1 = 20;          % First Stopband Attenuation (dB)
        Apass  = 1;           % Passband Ripple (dB)
        Astop2 = 20;          % Second Stopband Attenuation (dB)
        match  = 'stopband';  % Band to match exactly
        
        % Construct an FDESIGN object and call its BUTTER method.
        if Fstop2 == Fs/2
            h  = fdesign.highpass(Fstop1, Fpass1, Astop1, Apass, Fs);
            Hd = design(h, 'butter', 'MatchExactly', match);
        elseif Fstop1 == 0
            h  = fdesign.lowpass(Fpass2, Fstop2, Apass, Astop2, Fs);
            Hd = design(h, 'butter', 'MatchExactly', match);
        else
            h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                Astop2, Fs);
            Hd = design(h, 'butter', 'MatchExactly', match);
        end
        
        eval(['dataFilt', num2str(freqBandStop(band,1)), '_', num2str(freqBandStop(band,2)), '= zeros(size(dataFilt));'])
        
        for chan = 1:size(dataFilt,2)
            eval(['dataFilt', num2str(freqBandStop(band,1)), '_', num2str(freqBandStop(band,2)), '(:,chan)= filtfilt(Hd.sosMatrix, Hd.ScaleValues, dataFilt(:,chan));'])
            
            totalTime = toc(totalTimeTic);
            
            numChansComplete = chan + size(dataFilt,2)*(band-1);
            
            chanTime = totalTime/numChansComplete;
            timeLeft = chanTime * (size(dataFilt,2)-chan) + chanTime * size(dataFilt,2) * (size(freqBandPass,1)-band);
            
            clc;
            disp('*******************************')
            disp('Bandpass Filtering Data')
            disp([num2str(freqBandStop(band,1)), ' to ', num2str(freqBandStop(band,2)),' Hz Bandpass Filter. Chan: ', num2str(chan)])
            disp(['Time per chan: ', num2str(chanTime), ' sec'])
            disp(['Time spent: ', num2str(totalTime), ' sec'])
            disp(['Time Left: ', num2str(timeLeft), ' sec'])
            
        end % END FOR Low Pass Filter
        
        Header.Fband = [freqBandStop(band,:)];
        
        Header.bandParams.FiltType = 'butter';
        Header.bandParams.Fpass = freqBandPass;
        Header.bandParams.Fstop = freqBandStop;
        Header.bandParams.Astop1 = 20;
        Header.bandParams.Apass = 1;
        Header.bandParams.Astop2 = 20;
        Header.bandParams.match = 'stopband';
        
        
        band1 = num2str(freqBandStop(band,1));
        while length(band1)<3
            band1 = ['0', band1];
        end % END WHILE band1 name length
        
        band2 = num2str(freqBandStop(band,2));
        while length(band2)<3
            band2 = ['0', band2];
        end % END WHILE band1 name length
        
        bandOutputFullPath = [params.targetDir, fileName(1:end-4), '_', band1, '_', band2, '.mat'];
        
        Header.fileName = [fileName(1:end-4), '_', band1, '_', band2, '.mat'];
        Header.filePath = params.targetDir;
        
        eval(['data = dataFilt', num2str(freqBandStop(band,1)), '_', num2str(freqBandStop(band,2)), ';']);
        save(bandOutputFullPath, 'data', 'Header');
        clear data

        fileList{end+1} = Header.fileName;       
        
    end % END FOR num freq windows       
end % END IF params.bandpass

end % END FUNCTION

% EOF