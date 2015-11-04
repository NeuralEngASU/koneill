%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	nse2mat  (For Spike Events)
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
%		(SpikeData) is a structure containing the following fields
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
%			- exp2013-01-31_11_52_43_SE.mat     (yyyy-mm-dd_hh-mm-ss)
%				- inside .mat file there will be a structure named SpikeData
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = nse2mat( override )

	if nargin == 1
		fprintf('Override set to: %d\n', override);
	else
		override = 0;
		fprintf('Override set to: %d\n', override);
	end
	
	%***** Select Experimental Data to Load *****%
	[fileName, pathName] = uigetfile('*.nse', 'Select a Spike Event (.nse) file'); % Prompts the user to select *ONE* file.

	
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

	fileList = dir(fullfile(pathName, ['*', fileExt]));				% Gets the list of files with the selected file extension
	matExist = dir(fullfile(pathName, ['exp', expDate, '_SE.mat']));	% Gets the list of exp201x-xx-xx_xx-xx-xx_CSC.mat files

    % Reorders fileList
    tmpList = numel(fileList);
    
    for i = 1:tmpList
        fileStruct(i).name = ['SE', num2str(i), fileExt];
    end % END FOR
    
    fileList = fileStruct;
    
    
	if isempty(matExist) 
		% This line is is useless. It just checks to see if any errors will occur with reading in data.
        [~, ~, ~, ~, ~, ~] = Nlx2MatSpike( fullfile(pathName, [fileList(1).name]), [1,1,1,1,1], 1, 1, []);
    elseif override == 1
		fprintf('.mat file already exists. Overriding saved data.\n')
	else
		fprintf('.mat file already exists. If you want to replace it, set override = 1\nQuitting function.\n');
		return;
	end % END IF

	SpikeData.header 		= '';						% Header, saved as a string
	SpikeData.scNumbers 	= '';						% Unknown,
	SpikeData.cellNumbers 	= cell(numel(fileList), 1); % CellNumbers, somehow related to units. Use with spike sort 3D
	SpikeData.features 		= cell(numel(fileList), 1); % Contains 8 features: peak, valley
	SpikeData.samples 		= cell(numel(fileList), 1); % Contains the 32 samples for each spike (9 before/23 after, or similar)
	SpikeData.timeStamps 	= cell(numel(fileList), 1); % Contains the time stamps of when the spikes occurred
	SpikeData.expDate 		= expDate;					% Contains the date and time of the experiment
    SpikeData.ADBitVolt     = '';                       % Voltage Scaling unit

    % Gets the CSC time stamps in order to normalize the spike time stamps
    [timeBase] = Nlx2MatCSC( fullfile(pathName, 'CSC1.ncs'), [1,0,0,0,0], 0, 1, []);
    
	% Loops over the number of files in the directory
	for i = 1:numel(fileList)
	
		% Checks to see what file the loop is on. If the loop is NOT on the first file, complete the simple data parsing.
		% Else, start and initialize the CSCData structure.
		if i ~= 1
			[timeStamps, ~, ~, features, samples, ~] = Nlx2MatSpike( fullfile(pathName, [fileList(i).name]), [1,1,1,1,1], 1, 1, []);
			SpikeData.samples{i,1}    = samples;
			SpikeData.timeStamps{i,1} = timeStamps-timeBase(1);
			SpikeData.features{i,1}   = features;
		else
			[timeStamps, scNumbers, cellNumbers, features, samples, header] = Nlx2MatSpike( fullfile(pathName, [fileList(i).name]), [1,1,1,1,1], 1, 1, []);

			% Initializes the cell data values
			SpikeData.header = header;
			SpikeData.ScNumbers = scNumbers;
			SpikeData.cellNumbers = cellNumbers;

			
			SpikeData.samples{i,1} 	  = zeros(1, size(samples,2), size(samples,1) );
			SpikeData.timeStamps{i,1} = zeros(1, length(timeStamps));
			SpikeData.features{i,1}   = zeros(1, length(features));

			SpikeData.samples{i,1}    = samples;
			SpikeData.timeStamps{i,1} = timeStamps-timeBase(1);
			SpikeData.features{i,1}   = features;
		end % END IF
			
		fprintf('Channel %d/%d extracted.\n', i, numel(fileList))
			
	end % END FOR
    
    exprVolt = '([0-9]+$)';
    bitVolt = regexp(SpikeData.header{16}, exprVolt, 'tokens');	% Gets the file extension from the file name
    
    bitVolt = ['0.', bitVolt{1}{1}];
    
    SpikeData.ADBitVolt = str2double(bitVolt);
        
    
	save(fullfile(pathName,['exp', expDate, '_SE.mat']), 'SpikeData');	% Saves data in .mat file
	fprintf('\n********************\nCOMPLETED EXTRACTION\nFile output: %s\n\n', ['exp', expDate,'_SE.mat'])
end % END FUNCTION

% EOF