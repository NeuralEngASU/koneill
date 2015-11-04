%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SpikeSort
%   Desc: A short function for spike sorting the semester project data from
%         Artificial Neural Computation EEE Class.
%
%   Authors: Taylor Hearn, Kevin O'Neill, Denise Oswalt
%   Date: 2015.08.24
%
%   The purpose of this function is to prepare the data for the selected
%   sorting method.
%
%
%   Params:
%       Fs:         [Hz] The sampling rate of the data to be sorted.
%
%       numUnits:        (optional) The number of unique units in the data.  
%                        Or the number you want to classify. For new data, 
%                        there may be more or less than the ammount you
%                        specify. Specifying a number locks the sorter to 
%                        look for that many units.
%
%       method:    [1-3] Determins which method to use for sorting.
%                     (1): Use PCA + K-means, iterate k for best match or
%                          use defined numUnits. Randomly seed 20 times.
%                       2: ANN.
%
%       preProcess [0-3] Determine the ammount of preprocess to use before
%                        spike sorting.
%                       0: None. Use the raw data to spike sort.
%                       1: Filter. Highpass the signal. Use hpFs param to
%                          set the cuttoff frequency. Default is 250 Hz.
%                     (2): Remove line noise and filter. Uses detrend() and
%                          the Chronux toolbox's rmlinenoise() to remove
%                          offset and 60Hz signals. Afterwards filter.
%                       3: Remove line noise, auto-CAR, filter. Remove
%                          line noise/detrend, common average reference,
%                          filter.
%
%       hpFs        [Hz] The cuttoff frequency for a 4th order elliptic 
%                        high-pass filter.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = SpikeSort(Header, params, data)

%% Parse Input
if ~isfield(params, 'Fs');         Fs         = 500; else Fs         = params.Fs; end
if ~isfield(params, 'numSpikes');  numNeuron  =  -1; else numNeuron  = params.numNeuron; end
if ~isfield(params, 'method');     method     =   1; else method     = params.method; end
if ~isfield(params, 'preProcess'); preProcess =   2; else preProcess = params.preProcess; end
if ~isfield(params, 'hpFs');       hpFs       = 250; else hpFs       = params.hpFs; end


%% Preprocess

% Define a time vector. Currently in [seconds]
time = linspace(0, size(data,2)/Fs, size(data,2));
timeLimits = [time(1), time(end)];
timeLabel = 'Time, seconds';

voltageLimits = [-1.1*max(-1*data), 1.1*max(data)];
voltageLabel = ['Voltage, uV'];

switch(preProcess)
    case 0;
        % Do nothing
        anData = data;
        preProcessTitle = 'Raw data';
        
    case 1; % High-pass filter only
        % High-pass filter
        order = 10;   % Filter Order
        Fpass = hpFs; % First Passband Frequency [Hz]
        Astop = 60;   % First Stopband Attenuation [dB]
        Apass  = 1;   % Passband Ripple [dB]

        h = fdesign.highpass('N,Fp,Ast,Ap', order, Fpass, Astop, Apass);
        Hd = design(h, 'ellip');
        
        tmpData = filtfilt(Hd.sosMatrix,Hd.ScaleValues,tmpData);
        
        % Assign processed data
        anData = tmpData; 
        clear tmpData
        preProcessTitle = ['High-pass filter at ', num2str(hpFs), ' Hz'];
        
    case 2; % Detrend/remove line noise then high-pass filter
        % Detrend and remove 60Hz noise
        tmpData = detrend(data);
        
        paramsLine.tapers = [2,5];
        paramsLine.Fs     = Fs;
        paramsLine.fpass  = [0, Fs/2];
        paramsLine.pad    = 1;
        
        tmpData = rmlinesc(tmpData, paramsLine);
        
        % High-pass filter
        order = 10;   % Filter Order
        Fpass = hpFs; % First Passband Frequency [Hz]
        Astop = 60;   % First Stopband Attenuation [dB]
        Apass  = 1;   % Passband Ripple [dB]

        h = fdesign.highpass('N,Fp,Ast,Ap', order, Fpass, Astop, Apass);
        Hd = design(h, 'ellip');
        
        tmpData = filtfilt(Hd.sosMatrix,Hd.ScaleValues,tmpData);
        
        % Assign processed data
        anData = tmpData; 
        clear tmpData
        preProcessTitle = ['Detrended, High-pass filter at ', num2str(hpFs), ' Hz'];
        
    case 3; % Detrend/remove line noise, Auto-CAR?, high-pass filter
        % Detrend and remove 60Hz noise
        tmpData = detrend(data);
        
        paramsLine.tapers = [2,5];
        paramsLine.Fs     = Fs;
        paramsLine.fpass  = [0, Fs/2];
        paramsLine.pad    = 1;
        
        tmpData = rmlinesc(tmpData, paramsLine);
        
        % Auto-CAR/CMR
        % Look at data (or subsection) and determine the statistics of the noise
        % Generate a vector from noise using the same statistics
        % CAR data and the new noise vector.
        
        % High-pass filter
        order = 10;   % Filter Order
        Fpass = hpFs; % First Passband Frequency [Hz]
        Astop = 60;   % First Stopband Attenuation [dB]
        Apass  = 1;   % Passband Ripple [dB]

        h = fdesign.highpass('N,Fp,Ast,Ap', order, Fpass, Astop, Apass);
        Hd = design(h, 'ellip');
        
        tmpData = filtfilt(Hd.sosMatrix,Hd.ScaleValues,tmpData);
        
        % Assign processed data        
        anData = tmpData; 
        clear tmpData
        preProcessTitle = ['Detrended, Auto-CAR, High-pass filter at ', num2str(hpFs), ' Hz'];
    otherwise
end

% Plot sample of data

% subplot?
% sampleIdx?

plot(time, anData, 'k')
xlim(timeLimits)
ylim(voltageLimits)

xlabel(timeLabel, 'FontSize', 14)
ylabel(voltageLabel, 'FontSize', 14)
title(preProcessTitle, 'FontSize', 16)

figWidth  = 5;
figHeight = 3;
set(gcf, 'Position', [2, 2, figWidth, figHeight])



%% Find Spikes

% Look at Neural Engineering Code

if true
    
    sigma = 0:0.05:6;
    % find spikes above sigma and plot numSpikes vs sigma
    
    % surf (or whatever) plot the ISI for each sigma.
end

% Plot spikes. Look at SpikeWaveform.
% Use 'spaghetti' and 'patch-confidence' plots.
    

%% Sort Spikes

switch(method)
    case 0;
        % Plot PCA + cluster centers as well
            % Iterate for different # of clusters
                % Iterate for different seeds.
    case 1;
    case 2;
    otherwise
end

for u = 1:numUnits
    % Plot ISI for each Unit
    % Plot spike waveform for each unit
        % spaghetti and patch-confidence
end % END FOR

end % END FUNCTION

% EOF