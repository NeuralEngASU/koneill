%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	ncs2mat  (For CSCs)
%		Kevin O'Neill
%		Mario Capecchi Lab
%		20131102
%		v0.1
%		PI: Naveen Nagarajan
%		
%	Inputs:
%		override:
%			(0) Function will not save over existing experimental data (DEFAULT)
%			(1) Function will overwrite saved information. The information will be lost for all time.
%
%	Outputs:
%		There are no function outputs. But there is a saved file output.
%		(CSCData) is a structure containing the following fields
%			header:			Structure containing header information from channel 1
%			freq:			Integer containing the sampling frequency from channel 1
%			numValidSample:	Contains the number of valid samples(?) from channel 1
%			samples:		Contains a 2D matrix of samples. Rows = channels; Columns = voltages
%			timeStamps:		Contains a 2D matrix of timeStamps. Rows = channels; Columns = times
%
%	To Use:
%		Run function
%			- Function will quit if output file already exists
%			- To regenerate output file, set override = 1;
%		Data will be saved as:
%			- exp2013-01-31_11_52_43_CSC.mat     (yyyy-mm-dd_hh-mm-ss)
%				- inside .mat file there will be a structure named CSCData
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = ncs2mat( override )

	if nargin == 1
		fprintf('Override set to: %d\n', override);
	else
		override = 0;
		fprintf('Override set to: %d\n', override);
	end
	
	%***** Select Experimental Data to Load *****%
	[fileName, pathName] = uigetfile('*.ncs', 'Select a CSC (.ncs) file'); % Prompts the user to select *ONE* file.

	
	if isequal(fileName,0) % Cancels the script of no file was selected
		fprintf('User selected Cancel, program will end. DEBUG_Select_File\n');
		return;
	else
		exprExt  = '(\.[a-z]+$)';	% Expression to find the extension
		exprDate = '([0-9]+-[0-9]+-[0-9]+_[0-9]+-[0-9]+-[0-9]+)';	% Expression to find the date

		fileExt = regexp(fileName, exprExt, 'tokens');	% Gets the file extension from the file name
		expDate = regexp(pathName, exprDate, 'tokens');	% Gets the experiment date from the file path
		
        fileExt = fileExt{1}{1}; % Pulls the string outside of the cell
        expDate = expDate{1}{1};
        
		fprintf('%s files selected from: %s\nDate of experiment: %s\n', fileExt, pathName, expDate);
	end % END IF

    addpath(pathName);
    
	fileList = dir(fullfile(pathName, ['*', fileExt]));				% Gets the list of files with the selected file extension
	matExist = dir(fullfile(pathName, ['exp', expDate, '_CSC.mat']));	% Gets the list of exp201x-xx-xx_xx-xx-xx_CSC.mat files

    tmpList = numel(fileList);
    
    for i = 1:tmpList
        fileStruct(i).name = ['CSC', num2str(i), fileExt];
    end % END FOR
    
    fileList = fileStruct;
    
	if isempty(matExist)
		% This line is is useless. It just checks to see if any errors will occur with reading in data.
		[~, ~, ~, ~, ~, ~] = Nlx2MatCSC( fullfile(pathName, [fileList(1).name]), [1,1,1,1,1], 1, 1, 1);
    elseif override == 1
		fprintf('.mat file already exists. Overriding saved data.\n')
	else
		fprintf('.mat file already exists. If you want to replace it, set override = 1\nQuitting function.\n');
		return;
	end % END IF

	CSCData.header         = '';						% Header information, saved as a string	
	CSCData.freq           = '';						% Sampling rate, default = 32000 Hz
	CSCData.numValidSample = '';						% Unknown
	CSCData.samples        = cell(numel(fileList), 1);	% Samples, saved as a cell
	CSCData.timeStamps 	   = cell(numel(fileList), 1);	% Time stamps are for every 512th sample starting at 1. 
	CSCData.expDate        = expDate;					% Contains the date and time of the experiment
    CSCData.ADBitVolt      = [];                        % Contains the AD Bit Volt value
	
	% Loops over the number of files in the directory
	for i = 1:numel(fileList)
	
		% Checks to see what file the loop is on. If the loop is NOT on the first file, complete the simple data parsing.
		% Else, start and initialize the CSCData structure.
		if i ~= 1
			[timeStamps, samples] = Nlx2MatCSC( fullfile(pathName, [fileList(i).name]), [1,0,0,0,1], 0, 1, []);
			
            % Makes one timestamp for every sample
            TSFinal = [];
            for j = 1:length(timeStamps)
                
                tmpTime = linspace(timeStamps(j), timeStamps(j)+16000, 513);
                tmpTime = tmpTime(1:end-1);
                
                TSFinal = [TSFinal, tmpTime];
                
            end % END FOR
            
            TSFinal = TSFinal - TSFinal(1);
            
            CSCData.samples{i}    = samples(:)';
			CSCData.timeStamps{i} = TSFinal;
		else
			[timeStamps, channelNumbers, sampleFreq, numValidSample, samples, header] = Nlx2MatCSC( fullfile(pathName, [fileList(i).name]), [1,1,1,1,1], 1, 1, []);

			CSCData.header = header;
			CSCData.freq = sampleFreq;
			CSCData.channelNumbers = channelNumbers;
			
			% Initializes the cell data values
			CSCData.samples{i,1}    = zeros( 1, size(samples,2) * size(samples,1) );
            CSCData.timeStamps{i,:} = zeros( 1, length(timeStamps) * 512  );
            
            % Makes one timestamp for every sample
            TSFinal = [];
            for j = 1:length(timeStamps)
                
                tmpTime = linspace(timeStamps(j), timeStamps(j)+16000, 513);
                tmpTime = tmpTime(1:end-1);
                
                TSFinal = [TSFinal, tmpTime];
                
            end % END FOR
			
            TSFinal = TSFinal - TSFinal(1);
            
			CSCData.samples{i,1}    = samples(:)';
			CSCData.timeStamps{i,1} = TSFinal;
		end % END IF
			
		fprintf('Channel %d/%d extracted.\n', i, numel(fileList))
			
	end % END FOR
    
    exprVolt = '([0-9]+$)';
    bitVolt = regexp(CSCData.header{16}, exprVolt, 'tokens');	% Gets the file extension from the file name
    
    bitVolt = ['0.', bitVolt{1}{1}];
    
    CSCData.ADBitVolt = str2double(bitVolt);   
    
	
	save(fullfile(pathName,['exp', expDate, '_CSC.mat']), 'CSCData');	% Saves data in .mat file
	fprintf('\n********************\nCOMPLETED EXTRACTION\nFile output: %s\n\n', ['exp', expDate,'_CSC.mat'])
end % END FUNCTION

% EOF