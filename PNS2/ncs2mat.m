%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	ncs2mat  (For CSCs)
%		Kevin O'Neill
%		Mario Capecchi Lab
%		20140111
%		v1.0
%		PI: Naveen Nagarajan
%		
%	Inputs:
%		override:
%			(0) Function will not save over existing experimental data (DEFAULT)
%			(1) Function will overwrite saved information. The information
%			will be lost for all time (A very long time).
%
%	Outputs:
%		There are no function outputs. But there is a saved file output.
%		(CSCData) is a structure containing the following fields
%			header:			Structure containing header information from channel 1
%			freq:			Integer containing the sampling frequency from channel 1
%			numValidSample:	Contains the number of valid samples(?) from channel 1
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

function [] = ncs2mat( override, varargin )

	if nargin >= 1
		fprintf('Override set to: %d\n', override);
	else
		override = 0;
		fprintf('Override set to: %d\n', override);
    end % END IF
    
    dsFs = [];   % Downsample Frequency
    dsFlag = ''; % Downsample Flag
    
    if length(varargin) == 1
        dsFs = varargin{1};
        if ~(dsFs == 30000 || dsFs == 10000 || dsFs == 2000 || dsFs == 1000)
            fprintf('Need to choose a downsampled frequency: 30000, 10000, 2000, 1000\n')
            return;
        end % END IF
        fprintf('You must include the ''-ds'' flag id you wish to downsample.\nQuitting...\n')
        return;
    elseif length(varargin) == 2
        dsFs = varargin{1};
        if ~(dsFs == 30000 || dsFs == 10000 || dsFs == 2000 || dsFs == 1000)
            fprintf('Need to choose a downsampled frequency: 30000, 10000, 2000, 1000\n')
            return;
        end % END IF
        dsFlag = varargin{2};
        if ~strcmp(dsFlag, '-ds')
            fprintf('Downsample flag is incorrect, please use ''-ds''.\nQuitting...\n')
        else
            dsFlag = 1;
        end % END IF
    end % END IF
        
	
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
% 		[~, ~, ~, ~, ~, ~] = Nlx2MatCSC( fullfile(pathName, [fileList(1).name]), [1,1,1,1,1], 1, 1, 1);
    elseif override == 1
		fprintf('.mat file already exists. Overriding saved data.\n')
	else
		fprintf('.mat file already exists. If you want to replace it, set override = 1\nQuitting function.\n');
		return;
	end % END IF

	CSCData.header         = '';						% Header information, saved as a string	
	CSCData.freq           = '';						% Sampling rate, default = 32000 Hz
	CSCData.numValidSample = '';						% Unknown
	CSCData.timeStamps 	   = cell(numel(fileList), 1);	% Time stamps are for every 512th sample starting at 1. 
	CSCData.expDate        = expDate;					% Contains the date and time of the experiment
    CSCData.ADBitVolt      = [];                        % Contains the AD Bit Volt value
    CSCData.hdf5FileName   = [];                        % Contains the HDF5 file name
    CSCData.pathName       = pathName;                  % Contains the path informaton
    
    %% HDF5 File Name
    
    cscExt = '.h5';
    fsStr = num2str(dsFs);
    
    if dsFlag
        CSCMAT = ['exp', expDate, '_ds', fsStr, cscExt];
    else
        CSCMAT= ['exp', expDate, cscExt];
    end % END IF
    
    CSCData.hdf5FileName = CSCMAT;
	
    %% Extracts data and saves it to structure
    [timeStamps] = Nlx2MatCSC( fullfile(pathName, [fileList(i).name]), [1,0,0,0,0], 0, 1, []);
        
    timeStamps = timeStamps - timeStamps(1);
    CSCData.timeStamps = timeStamps; clear('timeStamps')
    
        
    [channelNumbers, sampleFreq, numValidSample, header] = Nlx2MatCSC( fullfile(pathName, [fileList(i).name]), [0,1,1,1,0], 1, 1, []);
    
    CSCData.header = header;                    clear('header')
    CSCData.freq = sampleFreq(1);               clear('sampleFreq')
    CSCData.channelNumbers = channelNumbers;    clear('channelNumbers')
    %%
    
    
    exprVolt = '([0-9]+$)';
    bitVolt = regexp(CSCData.header{16}, exprVolt, 'tokens');	% Gets the file extension from the file name
    
    bitVolt = ['0.', bitVolt{1}{1}];
    
    CSCData.ADBitVolt = str2double(bitVolt);   
    
	
	save(fullfile(pathName,['exp', expDate, '_CSC.mat']), 'CSCData');	% Saves data in .mat file
	fprintf('\n********************\nCOMPLETED HEADER EXTRACTION\nFile output: %s\n\n', ['exp', expDate,'_CSC.mat'])
    
%% Samples Extraction

cscExt = '.h5';
fsStr = num2str(dsFs);

if dsFlag
    CSCMAT = ['exp', expDate, '_ds', fsStr, cscExt];
else 
    CSCMAT= ['exp', expDate, cscExt];
end % END IF

fullCSCMAT = fullfile(pathName, CSCMAT);

% if exists(CSCMAT) && override
%     delete(CSCMAT)
% elseif exists(CSCMAT) && ~override
%     fprintf('.cscmat already exists. %s\nQuitting...\n', fullfile(pathName, CSCMAT))
%     return;
% end

% Creating filter object
if ~ isempty(dsFs)
    dsFs = CSCData.freq / round(CSCData.freq / dsFs);
    if CSCData.freq > dsFs
        N = 10;         % order
        Fpass = dsFs/4; % Passband frequency
        Apass = 1;      % Passband ripple (dB)
        Astop = 80;     % Stopband attenuation (dB)
        h = fdesign.lowpass('N,Fp,Ap,Ast',N,Fpass,Apass,Astop,CSCData.freq);
        Hd = design(h,'ellip');
    else
        fprintf('Cannot downsample to this frequency, %d.\nChoose a different value.\n', dsFs)
        return
    end % END IF
end % END IF

% Header information
Header.header         = CSCData.header;         % Text header from the .csc file
Header.freq           = 32000;                  % Sampling rate, default = 32000 Hz
Header.timeStamps 	  = CSCData.timeStamps;     % Time stamps are for every 512th sample starting at 1. 
Header.expDate        = expDate;				% Contains the date and time of the experiment
Header.ADBitVolt      = CSCData.ADBitVolt;      % Contains the AD Bit Volt value
Header.channelNumbers = CSCData.channelNumbers; % Unknown
Header.numValidSample = '';						% Unknown

clear('CSCData')

if ~isempty(dsFs)
    Header.freq = dsFs;
end % END IF

% save(fullCSCMAT, 'Header', '-v7.3')
% h5create(fullCSCMAT, 'Header', Header);

% Determining system memory to maximize data segments
[userview, ~] = memory;
maxBytes = userview.MaxPossibleArrayBytes;


% Calculating maximum channel samples to load into memory
segmentLength = maxBytes*0.25 / 8 / 512; % 8 bytes is the length of a double
% segmentLength = 512*1024/512;
dataLen = length(Header.timeStamps);

segNum = floor(dataLen/segmentLength) + 1;
segRemLen = rem(dataLen,segmentLength);

dataSize = dataLen*512*8*numel(fileList)/1024/1024/1024;

fprintf('\n********************\nSTARTING SAMPLE EXTRACTION\nFile output: %s\n\n', fullCSCMAT)
fprintf('NCS2MAT Parsing: %0.0f%% complete\nLoading %0.1f GB of data\n', 0, dataSize)
pause(5)


for i = 1:numel(fileList)
    
    tempChanName = ['CSC',num2str(i)];
    eval([tempChanName, ' = [];']);
%     save(fullCSCMAT,tempChanName,'-append')
    h5create(fullCSCMAT, ['/', tempChanName], [1, dataLen*512]); % Create a h5 file because of 32-bit windows restrictions
    
    for j = 1:segNum
        
        tempChanName = ['CSC',num2str(i)];
%         matObj = matfile(fullCSCMAT, 'Writable', true); matObj = matfile(fullCSCMAT);
        
        currDataSize = (dataLen*512*8*(i-1) + dataLen*512*8/segNum*(j)) / 1024/1024/1024;
        clc; fprintf('NCS2MAT Parsing: %0.0f%% complete\nLoading %0.1f GB of data\n', currDataSize/dataSize*100, dataSize)
        fprintf('Channel: %d\tSegment: %d\n',i,j);
                
        if j ~= segNum(end)
            [tempChan] = Nlx2MatCSC( fullfile(pathName, [fileList(i).name]), [0,0,0,0,1], 0, 2, [(j-1)*segmentLength + 1, (j-1)*segmentLength + segmentLength]);
        else
            [tempChan] = Nlx2MatCSC( fullfile(pathName, [fileList(i).name]), [0,0,0,0,1], 0, 2, [(j-1)*segmentLength + 1, (j-1)*segmentLength + segRemLen]);
        end % END IF
        
        [rows, cols] = size(tempChan);
        
        tempChan = reshape(tempChan', 1, rows*cols);
        tempChan = tempChan .* Header.ADBitVolt;       
        
        if ~isempty(dsFs)
%             eval([tempChan,' = filtfilt(Hd.sosMatrix,Hd.ScaleValues,double(',tempChan,'));'])
%             eval([tempChan,' = ',tempChan,'(1:',num2str(round(Header.freq/dsFs)),':end);'])
            tempChan = filtfilt(Hd.sosMatrix,Hd.ScaleValues,tempChan);
            tempChan = tempChan(1:round(Header.freq/dsFs):end);
            
        end
        
        eval([tempChanName, ' = tempChan; %clear(''tempChan'');']);
        
%         save(fullCSCMAT,tempChanName,'-append')
        h5write(fullCSCMAT, ['/', tempChanName], tempChan, [1, (j-1)*segmentLength + 1], [1, size(tempChan,2)]);
%         eval(['matObj.', tempChanName, ' = [matObj.', tempChanName, ',', tempChanName, '];']) 
        clear(tempChanName)
%         clear('matObj')
        
    end % END FOR
end % END FOR

fprintf('\n********************\nCOMPLETED SAMPLE EXTRACTION\nFile output: %s\n\n', fullCSCMAT)



end % END FUNCTION

% EOF