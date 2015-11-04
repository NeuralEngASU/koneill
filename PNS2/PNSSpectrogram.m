%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	MouseSpectrum
%		Kevin O'Neill
%		Mario Capecchi Lab
%		20140315
%		v0.1
%		PI: Naveen Nagarajan
%
%	Inputs:
%		CSCData: Data loaded from the expYYYY-MM-DD_hh-mm-ss_CSC.mat file
%											(nsc2mat output)
%
%       saveDate: Saves figures and mat files that corrospond to the
%                 calculated response.
%
%       channels: Which channels 
%
%       calcOption: [1] Default. Calculates the FFT over the entire given 
%                   event and timeBounds arguments.
%
%                   [2] Calculates the FFT using a series of 400ms windows
%                   which occur every 50ms. Data will be truncated if the 
%                   timeBounds given cannot be divided 50ms. TimeBounds 
%                   must also be greater than 400ms
%
%	Outputs:
%		chanFFT:    A [numChan x NFFT x events] matrix that contains the
%                   FFT calculation.
%
%       Freq:       A vector containing the frequency over which the FFT
%                   calculated.
%
%	To Use:
%		Run function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ chanSpec, freq ] = MouseSpectrogram( CSCData, saveData, channels, timeBounds, events)

% Ensures that there is at least one arbitrary event
if ~isempty(events)
    numEvents = size(events);
else
    numEvents = 1;
    events = 0;
end % END IF

% Ensures that some channels are selected
if isempty(channels)
    channels = [1:16];
    numChans = length(channels);
else
    numChans = length(channels);
end % END IF

% Initializes output variables
chanSpec = [];
freq = [];

% Loops over the number of events given. If no event was given, then loop
% over the single time given.
for i = 1:numEvents
    startTime = events(1) + timeBounds(1);	% Seconds
    endTime   = events(1) + timeBounds(2);	% Seconds
    
    startTime = startTime * 1e6;			% micro seconds
    endTime   = endTime   * 1e6;           % micro seconds
    
    timeStamps = CSCData.timeStamps; % vector containing the timestamps for all channels
    timeIdx = timeStamps > startTime & timeStamps < endTime;
    
    L = length(timeIdx); % Length of the sample being analyzed
    NFFT = nextpow2(L);  % Next power of 2 greater than L
    
    for chan = 1:numChans
        
        tempChanName = ['CSC',num2str(channels(chan))];
        chanData(i,:) = h5read(CSCData.hdf5FileName, ['/', tempChanName], timeIdx(1), timeIdx(end) - timeIdx(1));
        
    end % END FOR
    
        
    movingWin = [];
    
    %                                              [TW, 2TW-1] % W = bandwidth, T = time duration
    %                          [Time-bandwidth product, number of tapers]
    params.tapers = [CSCData.freq/2 * diff(timeBounds), 2*CSCData.freq/2 * diff(timeBounds)-1]; 
    params.pad = 0; % Default padding
    params.Fs = CSCData.freq;
    params.fpass = [minFreq, maxFreq]; % unknown if the minimum and maximum pass frequency should be controlled by the user
                                       % needs testing in order to find
                                       % appropriot values
    
    % Please look at the chronux manual for the mtspecrtrogram
    % Might need a FOR loop for the chanData variable.
    [S, t, f, Serr] = mtspecgramc( chanData, movingWin, params);

    
end % END FOR


%%% Plot Here %%%
% Most likley will need a FOR loop

    shortName = ['SpectGram_', SpikeData.expDate, '_chan', num2str(k)];
    
    if saveData
        print('-dpng', ['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.png'], '-r100');
        saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.fig']);
        %			saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.png']);
%         saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.eps'],'epsc2');
        saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.pdf']);

    end % END IF

end % END FUNCTION

% EOF