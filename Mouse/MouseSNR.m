%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	MouseSNR  
%		Kevin O'Neill
%		Mario Capecchi Lab
%		20131103
%		v0.1
%		PI: Naveen Nagarajan
%		
%	Inputs:
%		CSCData: Data loaded from the expYYYY-MM-DD_hh-mm-ss_CSC.mat file
%											(nsc2mat output)
%       SpikeData
%
%	Outputs:
%		SNR = pk2pk / 2*stdev(noise); 
%
%	To Use:
%		Run function
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sigNoise] = MouseSNR( CSCData, SpikeData, channel, noiseTime )
	
	if nargin < 3
		fprintf('Not enough input arguments to MouseSNR. Quitting program. DEBUG_nargin\n');
		return
	end %END IF
    
    noiseTime = noisetime * 1e6; % Covnerts seconds to microseconds
    
    % Loop over given hannels
    for i = 1:length(channel)
        
        % Select channel and pre-allocate memory
        chan = channel(i);
        SNRTemp = zeros(1, numel(SpikeData.samples{chan}));
        
        % Loop over the number of spikes within a channel
        for j = 1:numel(SpikeData.samples{chan})
            
            % Find the peak-to-peak voltage of a spike
            voltData = SpikeData.samples{chan,j};
            pk2pk = max(voltData) - min(voltData);  % There may be a better way to find pk2pk as this max()-min() feels inelegant
            
            % Find times for baseline
            timeIdx = find((CSCData.timeStamps >= noiseTime(1) & CSCData.timeStamps <= noiceTime(2)) == 1);
            
            % Load baseline data
            tempChanName = ['CSC',num2str(chan)];
            noiseData = h5read(CSCData.hdf5FileName, ['/', tempChanName], timeIdx(1), length(timeIdx)); 
            
            % Calculate SNR
            SNRTemp(j) = pk2pk / 2*stdev(noiseData);
            
        end % END FOR
        
        % Average and Standard Error SNR for each channel
        sigNoise(i,1) = mean(SNRTemp);
        sigNoise(i,2) = stdev(SNRTemp) / sqrt(length(SNRTemp));
        
    end % END FOR
	
end % END FUNCTION
% EOF