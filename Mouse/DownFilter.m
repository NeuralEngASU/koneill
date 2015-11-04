%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	DownFilter 
%		Kevin O'Neill
%		Mario Capecchi Lab
%		20140411
%		v1.0
%		PI: Naveen Nagarajan
%
%	Inputs:
%		override:
%			(0) Function will not save over existing experimental data (DEFAULT)
%			(1) Function will overwrite saved information. The information
%			will be lost for all time (A very long time).
%
%       CSCData:
%           Data loaded from the expYYYY-MM-DD_hh-mm-ss_CSC.mat file
%											(nsc2mat output)
%
%       dsFlag:
%           Boolean. Determines whether the samples are downsampled
%
%       dsFs:
%           Integer. The sampling frequency for the downsampled data. The
%           maximum frequency after downsampling is 1/2 dsFs.
%
%       flitFlag
%           Boolean. Determines whether the downsampled data is then
%           filtered. Or if dsFlag is false, filters the given data.
%
%       flitType:
%           String. 'highpass' or 'lowpass' are the only inputs. These
%           allow the user to control the filtering object.
%
%           In order to show spikes use a 'highpass' filter with a value of
%           500 Hz.
%
%           For LFP data, please use a 'lowpass' with a value of 500 Hz. 
%
%       filtFs:
%           Double/integer. The filter frequency used in the filters.
%
%	Outputs:
%		There are no function outputs. But there is a saved file output.
%		(CSCData) is a structure containing the following fields
%			header:			Structure containing header information from channel 1
%			freq:			Integer containing the sampling frequency from channel 1
%			numValidSample:	Contains the number of valid samples(?) from channel 1
%			timeStamps:		Contains a 2D matrix of timeStamps. Rows = channels; Columns = times
%
%       HDF5 File:
%           Contains the saved samples data after they have been
%           downsampled and/or filtered. 
%
%	To Use:
%		Run function
%			- Function will quit if output file already exists
%			- To regenerate output file, set override = 1;
%		CSCData will be saved as:
%			- exp2013-01-31_11_52_43_CSC.mat     (yyyy-mm-dd_hh-mm-ss)
%				- inside .mat file there will be a structure named CSCData
%		HDF5 data will be saved as:
%			- exp2013-01-31_11_52_43_CSC_dsFREQ_filtFREQ.h5     (yyyy-mm-dd_hh-mm-ss)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = DownFilter( override, CSCData, dsFlag, dsFs, filtFlag, filtType, filtFs)

if dsFlag && ~(dsFs == 30000 || dsFs == 10000 || dsFs == 2000 || dsFs == 1000)
    fprintf('Need to choose a downsampled frequency: 30000, 10000, 2000, 1000\n')
    return;
end % END IF

if filtFlag && ~(isequal(filtType, 'highpass') || isequal(filtType, 'lowpass'))
    fprintf('Need to select a filtered type: highpass or lowpass.\n')
    return;
end % END IF

if filtFs > 0.5*dsFs
    fprintf('Filtered frequency needs to be less than or equal to 1/2 * samples (dsFs).\n')
    return;
end % END IF


%% HDF5 File Name

cscExt = '.h5';
dsStr = num2str(dsFs);
fsStr = num2str(filtFs);

CSCMAT = ['exp', expDate];

if dsFlag
    CSCMAT = ['exp', expDate, '_ds', dsStr];
end

if filtFlag
    CSCMAT = [CSCMAT, '_', filtType, fsStr];
end

CSCData.hdf5FileName = [CSCMAT, cscExt];

%% New CSCData Save File

fullCSCh5 = fullfile(CSCData.pathName, [CSCMAT, cscExt]);
fullCSCMAT = fullfile(CSCData.pathName, [CSCMAT, 'mat']);

Header.freq = CSCData.freq; % Temporary variable for the frequency of the output.

%% Creating downsample object
if dsFlag
    dsFs = CSCData.freq / round(CSCData.freq / dsFs);
    if CSCData.freq > dsFs
        N = 10;         % order
        Fpass = dsFs/4; % Passband frequency
        Apass = 1;      % Passband ripple (dB)
        Astop = 80;     % Stopband attenuation (dB)
        h = fdesign.lowpass('N,Fp,Ap,Ast',N,Fpass,Apass,Astop,CSCData.freq);
        Hd = design(h,'ellip');
        
        Header.freq = dsFs;
    else
        fprintf('Cannot downsample to this frequency, %d.\nChoose a different value.\n', dsFs)
        return
    end % END IF
end % END IF


%% Creating filter object
if filtFlag
    if isequal(filtType, 'highpass')
        N = 10;         % order
        Fpass = filtFs; % Passband frequency
        Apass = 1;      % Passband ripple (dB)
        Astop = 80;     % Stopband attenuation (dB)
        h = fdesign.highpass('N,Fp,Ap,Ast',N,Fpass,Apass,Astop,Header.freq);
        Hf = design(h,'ellip');
    elseif isequal(filtType, 'lowpass')
        N = 10;         % order
        Fpass = filtFs; % Passband frequency
        Apass = 1;      % Passband ripple (dB)
        Astop = 80;     % Stopband attenuation (dB)
        h = fdesign.lowpass('N,Fp,Ap,Ast',N,Fpass,Apass,Astop,Header.freq);
        Hf = design(h,'ellip');        
        
    else
        fprintf('There was an error _DEBUG_.\n')
        return
    end % END IF
end % END IF


%% Find usable memory size
% Determining system memory to maximize data segments
[userview, ~] = memory;
maxBytes = userview.MaxPossibleArrayBytes;


% Calculating maximum channel samples to load into memory
segmentLength = maxBytes*0.4 / 8 / 512; % 8 bytes is the length of a double
% segmentLength = 512*1024/512;
dataLen = length(Header.timeStamps);

segNum = floor(dataLen/segmentLength) + 1;
segRemLen = rem(dataLen,segmentLength);

dataSize = dataLen*512*8*16/1024/1024/1024;

%% Downsample and Filter
fprintf('\n********************\nSTARTING SAMPLE EXTRACTION\nFile output: %s\n\n', fullCSCMAT)
fprintf('DOWNSAMPLE Parsing: %0.0f%% complete\nLoading %0.1f GB of data\n', 0, dataSize)
pause(5)

if dsFlag
    CSCData.timeStamps = CSCData.timeStamps(1:round(CSCData.freq/dsFs):end);
end % END IF

save(fullCSCMAT, 'CSCData');	% Saves data in .mat file

maxChans = 16; % This may need to be changed to be non-hardcoded

for i = 1:maxChans   
    
    tempChanName = ['CSC',num2str(i)];
    eval([tempChanName, ' = [];']);

    % Creates new hdf5 file
    h5create(fullCSCh5, ['/', tempChanName], [maxChans, dataLen*512]); % Create a h5 file because of 32-bit windows restrictions

    
    for j = 1:segNum
        
        tempChanName = ['CSC',num2str(i)];
        %         matObj = matfile(fullCSCMAT, 'Writable', true); matObj = matfile(fullCSCMAT);
        
        currDataSize = (dataLen*512*8*(i-1) + dataLen*512*8/segNum*(j)) / 1024/1024/1024;
        clc; fprintf('DownFilter Parsing: %0.0f%% complete\nLoading %0.1f GB of data\n', currDataSize/dataSize*100, dataSize)
        fprintf('Channel: %d\tSegment: %d\n',i,j);
        
        if j ~= segNum(end)
            samples = [(j-1)*segmentLength + 1, (j-1)*segmentLength + segmentLength];
            [tempChan] = h5read(fullfile(CSCData.pathName, CSCData.hdf5FileName), ['/', tempChanName], samples(1), samples(1)-samples(2));
        else
            samples = [(j-1)*segmentLength + 1, (j-1)*segmentLength + segRemLen];
            [tempChan] = h5read(fullfile(CSCData.pathName, CSCData.hdf5FileName), ['/', tempChanName], samples(1), samples(1)-samples(2));
        end % END IF
        
        if dsFlag
            tempChan = filtfilt(Hd.sosMatrix,Hd.ScaleValues,tempChan);
            tempChan = tempChan(1:round(CSCData.freq/dsFs):end);
        end
        
        if filtFlag
            tempChan = filtfilt(Hf.sosMatrix,Hf.ScaleValues,tempChan);
        end
                
        h5write(fullCSCMAT, ['/', tempChanName], tempChan, [1, (j-1)*segmentLength + 1], [1, size(tempChan,2)]);
        clear(tempChanName)
        clear('tempChan')

        
    end % END FOR
end % END FOR

fprintf('\n********************\nCOMPLETED DOWNSAMPLE AND FILTERING\nFile output: %s\n\n', fullCSCMAT)



end % END FUNCTION

% EOF
