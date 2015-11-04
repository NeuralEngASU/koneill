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
%       calcOption: [1] Default. Calculates the FFT over the entire given 
%                   event and timeBounds arguments.
%
%                   [2] Calculates the FFT using a series of 400ms windows
%                   which occur over a 50ms trains. Data will be truncated
%                   if the timeBounds given cannot be divided 50ms.
%                   TimeBounds must also be greater than 400ms
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

function [ chanFFT, freq ] = MouseSpectrum( CSCData, saveData, channels, timeBounds, events, nameMod)


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
chanFFT = [];
freq = [];

% Loops over the number of events given. If no event was given, then loop
% over the single time given.
for i = 1:numEvents
    startTime = events(1) + timeBounds(1);	% Seconds
    endTime   = events(1) + timeBounds(2);	% Seconds
    
    startTime = startTime * 10e5;			% micro seconds
    endTime   = endTime   * 10e5;           % micro seconds
    
    timeStamps = CSCData.timeStamps; % vector containing the timestamps for all channels
    timeIdx = timeStamps > startTime & timeStamps < endTime;
    
    L = length(timeIdx); % Length of the sample being analyzed
    NFFT = 2^nextpow2(L);  % Next power of 2 greater than L
    
    for chan = 1:numChans
        
        tempChanName = ['CSC',num2str(channels(chan))];
        chanData(i,:) = h5read(CSCData.hdf5FileName, ['/', tempChanName], timeIdx(1), timeIdx(end) - timeIdx(1));
        
    end % END FOR
    
    % Please look at Chronux mtspectrumc
    chanFFT(:,:,i) = fft(chanData, NFFT, 1);
    
end % END FOR

freq = CSCData.freq/2 * linspace(0, 1, NFFT/2+1);

plot(f,2*abs(Y(1:NFFT/2+1)))

if saveData
    for i = 1:size(chanFFT, 2)
        for j = 1:size(chanFFT, 3)
            
            % Plot
            
            
            
        end % END IF
    end % END IF
end % END IF
end % END FUNCTION

% EOF