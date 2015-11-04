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
%       filtFlag
%           Boolean. Determines whether the downsampled data is then
%           filtered. Or if dsFlag is false, filters the given data.
%
%       filtType:
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

function [] = PNSDownFilter( override, nsxFilePath, varargin )

% Defaults
DSFLAG = 0;
DSFS = 2000;
FILTFLAG = 0;
FILTTYPE = 'lowpass';
FILTFS = 1000;
TIMESTART = 0; % Seconds
TIMESPAN = 100; % Seconds

% parse varargin
for ii = 1:2:length(varargin)
    if ~exist(upper(varargin{ii}), 'var')
        fprintf('Unknown option entry: %s\n', varargin{ii})
        return;
    else
        eval([upper(varargin{ii}) '=varargin{ii+1};']);
    end
end
dsFlag = DSFLAG;
dsFs = DSFS;
filtFlag = FILTFLAG;
filtType = FILTTYPE;
filtFs = FILTFS;
timeStart = TIMESTART;
timeSpan = TIMESPAN;

dsFsStr = num2str(dsFs);
filtFsStr = num2str(filtFs);

clear DSFLAG DSFS FILTFLAG FILTTYPE FILTFS

if dsFlag && ~(dsFs == 30000 || dsFs == 10000 || dsFs == 2000 || dsFs == 1000)
    fprintf('Need to choose a downsampled frequency: 30000, 10000, 2000, 1000\n')
    return;
end % END IF

if filtFlag && ~(isequal(filtType, 'highpass') || isequal(filtType, 'lowpass'))
    fprintf('Need to select a filtered type: highpass or lowpass.\n')
    return;
end % END IF

if filtFs > 0.5*dsFs
    fprintf('Filtered frequency needs to be less than 1/2 * samples (dsFs).\n')
    return;
end % END IF

if ~dsFlag && ~filtFlag
    fprintf('Need to specify dsFlag or filtFlag.\n');
    return;
end % END IF

%% HDF5 File Name

pathExpr = '(.+)?\\';
pathName = regexp(nsxFilePath, pathExpr, 'Tokens');
pathName = pathName{1}{1};

nameExpr = '.+\\(.+\..+$)';
fileName = regexp(nsxFilePath, nameExpr, 'Tokens');
fileName = fileName{1}{1};

fileNameMat = fileName(1:end-4);

if dsFlag
    fileNameMat = [fileNameMat, '_ds', dsFsStr];
end

if filtFlag
    fileNameMat = [fileNameMat, '_', filtType, filtFsStr];
end

fileNameMat = [fileNameMat, fileName(end-3:end), 'mat'];

%% New Header Info

Header = openNSx(nsxFilePath);
Header = Header.MetaTags;

%% Creating downsample object
if dsFlag
    dsFs = Header.SamplingFreq / round(Header.SamplingFreq / dsFs);
    if Header.SamplingFreq > dsFs
        N = 10;         % order
        Fpass = dsFs/2; % Passband frequency
        Apass = 1;      % Passband ripple (dB)
        Astop = 80;     % Stopband attenuation (dB)
        h = fdesign.lowpass('N,Fp,Ap,Ast',N,Fpass,Apass,Astop,Header.SamplingFreq);
%         h = fdesign.lowpass('N,F3dB',N, Fpass,Header.SamplingFreq);
        Hd = design(h,'ellip');
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
%         h = fdesign.highpass('N,Fp,Ap,Ast',N,Fpass,Apass,Astop,Header.SamplingFreq);
        h = fdesign.highpass('N,F3dB',N, Fpass,Header.SamplingFreq);
        Hf = design(h,'butter');
    elseif isequal(filtType, 'lowpass')
        N = 10;         % order
        Fpass = filtFs; % Passband frequency
        Apass = 1;      % Passband ripple (dB)
        Astop = 80;     % Stopband attenuation (dB)
%         h = fdesign.lowpass('N,Fp,Ap,Ast',N,Fpass,Apass,Astop,Header.SamplingFreq);
        h = fdesign.lowpass('N,F3dB',N, Fpass,Header.SamplingFreq);
        Hf = design(h,'butter');        
        
    elseif isequal(filtType, 'bandpass')
        N = 10;         % order
        Flowpass = filtFs(1); % Lower passband frequency
        Fhighpass = filtFs(2); % Higher passband frequency
        Apass = 1;      % Passband ripple (dB)
        Astop = 80;     % Stopband attenuation (dB)
%         h = fdesign.bandpass('N,Fp1,Fp2,Ast1,Ap,Ast2',N,Flowpass,Fhighpass,Astop,Apass,Astop,Header.SamplingFreq);
%         h = fdesign.lowpass('N,F3dB',N, Fpass,Header.SamplingFreq);
%         Hf = design(h,'ellip');      
    else
        fprintf('There was an error _DEBUG_.\n')
        return
    end % END IF
end % END IF


%% Find usable memory size
% % Determining system memory to maximize data segments
% [userview, ~] = memory;
% maxBytes = userview.MaxPossibleArrayBytes;
% 
% 
% % Calculating maximum channel samples to load into memory
% segmentLength = maxBytes*0.4 / 8 / 512; % 8 bytes is the length of a double
% % segmentLength = 512*1024/512;
% dataLen = length(Header.timeStamps);
% 
% segNum = floor(dataLen/segmentLength) + 1;
% segRemLen = rem(dataLen,segmentLength);
% 
% dataSize = dataLen*512*8*16/1024/1024/1024;

%% Downsample and Filter
sizeFile = dir(nsxFilePath);
sizeFile = sizeFile.bytes;

sizeChan = sizeFile/Header.ChannelCount;

fprintf('\n********************\nSTARTING SAMPLE EXTRACTION\nFile output: %s\n\n', fullfile(pathName,fileNameMat))
fprintf('DOWNFILTER Parsing: %0.0f%% complete\nLoading %0.1f GB of data\n', 0, sizeFile/(1024^3))
pause(1)

save(fullfile(pathName, fileNameMat), 'Header')

maxChans = 96; %Header.ChannelCount; % This may need to be changed to be non-hardcoded
if(timeStart == 0)
    timeStartSamp = 1; 
else
    timeStartSamp = ceil(Header.SamplingFreq*timeStart); 
end

if(timeSpan == -1)
    timeSpanSamp = Header.Samples; 
else
    timeSpanSamp = ceil(timeSpan * Header.SamplingFreq + timeStartSamp); 
end

disp('Loading Data')
loadTime = tic;
nsxData = openNSx(nsxFilePath, 'read', ['c:1:', num2str(maxChans)], ['t:', num2str(timeStartSamp), ':', num2str(timeSpanSamp)];
fprintf('Data Loaded: %f seconds\n', toc(loadTime))
timeSpent = tic;
for ii = 1:maxChans   
    
%     nsxData = openNSx(nsxFilePath, 'read', ['c:' num2str(ii)]);
        
    tempChan = double(nsxData.Data(ii,:));
    
    if filtFlag
        filtTime = tic;
        disp('Filter Start')
        tempChan = filtfilt(Hf.sosMatrix,Hf.ScaleValues,tempChan);
        filtTimeOut = toc(filtTime);
        fprintf('Filter End: %f seconds\n', filtTimeOut)
    end
    
    if dsFlag
        dsTime = tic;
        disp('DownSample Start')
        tempChan = filtfilt(Hd.sosMatrix,Hd.ScaleValues,tempChan);
        tempChan = tempChan(1:round(Header.SamplingFreq/dsFs):end);
        dsTimeOut = toc(dsTime);
        fprintf('DownSample End: %f seconds\n', dsTimeOut)
    end
    
    tempChanName = ['C', num2str(ii)];
    
    eval([tempChanName,' = tempChan;']);
    save(fullfile(pathName, fileNameMat), tempChanName, '-append')
    clear('tempChan')

    clc; fprintf('DownFilter Parsing: %0.0f%% complete\nLoading %0.1f GB of data\n', (sizeChan*ii)/sizeFile*100, sizeFile/(1024^3))
    fprintf('Channel: %d\n',ii);    
    fprintf('Time spent: %f seconds\n', toc(timeSpent))
end % END FOR

fprintf('\n********************\nCOMPLETED DOWNSAMPLE AND FILTERING\nFile output: %s\n\n',  fullfile(pathName,fileNameMat))



end % END FUNCTION

% EOF