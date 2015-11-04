%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% blocksev2mat.m  (based on SEV2mat.m from TDT)
%   Author: Kevin O'Neill
%   Date: 2015/07/29
%   Desc:
%       Used to extract data from the recorded .sev files. Each .sev file
%       is assigned to one channel and contains a 40 byte header followed
%       by data.
%
%       sourceDir: the source folder for the .sev files
%
%       targetDir: the target folder for the .mat file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Header = blocksev2mat(Header, params)
%% Initialize variables

sourceDir = Header.info.blockPath; % The source directory of the *.sev files (block)
targetDir = ''; % The location where you want to save the extracted data

if ~isfield(params, 'targetDir'); targetDir = uigetdir('', 'Select Output Folder'); else targetDir = params.targetDir; end 

startTime = '';

%% Load tank/block time information

if isempty(startTime)
    tsqList =  dir(fullfile(sourceDir, '*.tsq'));

    if isempty(tsqList)
        error('Cannot find block data in source directory.');
    end % END IF isempty(tsqList)
    
    % Open Block data file
    tsq = fopen(fullfile(sourceDir, tsqList(1).name), 'r');
    
    % allocate variables
    tsqStartTimeStamp = 0;
    tsqStopTimeStamp = 0;
    discard = 0;
    count = 0;
    
    % Find the first initialized startTime
    while (tsqStartTimeStamp == 0 || count >=100)
        discard = fread(tsq, 1, 'long');
        discard = fread(tsq, 1, 'long');
        discard = fread(tsq, 1, 'long');
        discard = fread(tsq, 1, 'uint16');
        discard = fread(tsq, 1, 'uint16');
        tsqStartTimeStamp = fread(tsq, 1, 'double');
        discard = fread(tsq, 1, 'uint64');
        discard = fread(tsq, 1, 'long');
        discard = fread(tsq, 1, 'float');
        
        count = count+1;
    end % END WHILE tsqTimeStamp
    
    % Go to the start of the last timestamp
    fseek(tsq, -24, 'eof');
    
    % Extract the last timestamp
    tsqStopTimeStamp = fread(tsq, 1, 'double');
        
    fclose(tsq);
   
    timeRef = datenum('1970', 'yyyy'); % Setup a reference date (Jan-1 1970)
    startTimeMatlab = timeRef + tsqStartTimeStamp / 8.64e4; % Convert the tsqStartTimeStamp into days and add the days to the reference date
    startTimeMatlab = startTimeMatlab - 7/24; % Arizona time is GMT-7:00
    startTimeMatlabString = datestr(startTimeMatlab, 'yyyymmdd HH:MM:SS'); % Convert to string
    
    stopTimeMatlab = timeRef + tsqStopTimeStamp / 8.64e4; % Convert the tsqStopTimeStamp into days and add the days to the reference date
    stopTimeMatlab = stopTimeMatlab - 7/24; % Arizona time is GMT-7:00
    stopTimeMatlabString = datestr(stopTimeMatlab, 'yyyymmdd HH:MM:SS'); % Convert to string
    
    duration = stopTimeMatlab - startTimeMatlab;
    durationString = datestr(duration, 'HH:MM:SS');
    
    Header.info.dateStart = startTimeMatlabString(1:8);
    Header.info.dateStop = stopTimeMatlabString(1:8);
    Header.info.startTime = startTimeMatlabString(10:end);
    Header.info.stopTime = stopTimeMatlabString(10:end);
    Header.info.duration = durationString;
else
    Header.info.dateStart = '';
    Header.info.dateStop = '';
    Header.info.startTime = startTime;
    Header.info.stopTime = -1;
    Header.info.duration = -1;
end % END IF

Header.info.targetDir = targetDir;

startTime = Header.info.startTime;

%% Extract Header and Test for the length of the sev files

fileList = dir(fullfile(sourceDir, '*.sev'));

if length(fileList) < 1
    warning(['No .sev files found in: ', sourceDir])
end % END IF length(fileList) < 1

% List of allowed formats
allowedFormats = {'single','int32','int16','int8','double','int64'};
allowedFormatsSize = [16, 32, 16, 8, 32, 64];
% Output the Header once to a .mat file
headerCount = 0;

% Extract the filename for the data.
exprStr = '([A-Za-z0-9\_]+)\_DSP[0-9]\_Ch[0-9]\.[sev]';
fileName = regexp(fileList(1).name, exprStr, 'Tokens');
fileName = fileName{1}{1};

outputName = [fileName, '.mat'];
outputPath = fullfile(targetDir, outputName);

filePath = fullfile(sourceDir, fileList(1).name);
FID = fopen(filePath, 'rb');

if FID < 0
    warning(['Cannot open: ' filePath])
    return
end % END IF FID < 0

tmpHeader = struct();

tmpHeader.fileSizeBytes = fread(FID,1,'uint64');
tmpHeader.fileType      = char(fread(FID,3,'char')');
tmpHeader.fileVersion   = fread(FID,1,'char');

% Extrct Header information
if tmpHeader.fileVersion < 3
    % Event name
    if tmpHeader.fileVersion == 2
        tmpHeader.eventName  = char(fread(FID,4,'char')');
    else
        tmpHeader.eventName  = fliplr(char(fread(FID,4,'char')'));
    end % END IF fileVersion
    
    tmpHeader.channel     = fread(FID, 1, 'uint16'); % Current channel
    tmpHeader.numChan     = fread(FID, 1, 'uint16'); % Number of channels
    tmpHeader.numByteSamp = fread(FID, 1, 'uint16'); % Number of bytes per sample
    reserved              = fread(FID, 1, 'uint16'); % Reserved segment
    
    % Data format in lower four bits
    tmpHeader.dataFormat = allowedFormats{bitand(fread(FID, 1, 'uint8'),7)+1};
    
    % Items used to calculate Fs
    tmpHeader.decimate   = fread(FID, 1, 'uint8');
    tmpHeader.rate       = fread(FID, 1, 'uint16');
    
    % Reserved tags
    reserved = fread(FID, 1, 'uint64');
    reserved = fread(FID, 2, 'uint16');
    
else
    error(['Unknown file version: ' num2str(tmpHeader.fileVersion)]);
end % END IF fileVersion

fclose(FID);

% Determine sampling rate
if tmpHeader.fileVersion > 0
    tmpHeader.Fs = 2^(tmpHeader.rate)*25000000/2^12/tmpHeader.decimate;
    %         exists = isfield(tmpData, streamHeader.eventName);
else
    tmpHeader.dForm = 'single';
    tmpHeader.Fs = 0;
    s = regexp(file_list(ii).name, '_', 'split');
    tmpHeader.eventName = s{end-1};
    tmpHeader.channelNum = str2double(regexp(s{end}, '\d+', 'match'));
    warning('%s has empty header; assuming %s ch %d format %s and fs = %.2f\nupgrade to OpenEx v2.18 or above\n', ...
        file_list(ii).name, tmpHeader.eventName, ...
        tmpHeader.channelNum, tmpHeader.dataFormat, 24414.0625);
    
    exists = 1;
    %data.(tmpHeader.eventName).fs = tmpHeader.fs;
    tmpHeader.Fs = 24414.0625;
end % END IF fileVersion > 0

% numSamp = tmpHeader.fileSizeBytes/ tmpHeader.numByteSamp - 10;
% numSec = numSamp / tmpHeader.Fs;

% if numSec <= time2Extract(2)
%     error('The time-of-interest range is beyond the scope of this file');
% end % END IF

Header.info.sevFs = tmpHeader.Fs;
Header.info.fileVersion = tmpHeader.fileVersion;
% Header.info.numSamp = numSamp;
% Header.numSampSegment = diff(time2Extract)*Header.Fs;
Header.info.dataFormat = 'double';
Header.info.outputPath = outputPath;
Header.info.outputName = outputName;

%% Extract data
for ii = 1:length(fileList)
    
    filePath = fullfile(sourceDir, fileList(ii).name);
    FID = fopen(filePath, 'rb');
    
    if FID < 0
        warning(['Cannot open: ' filePath])
        return
    end % END IF FID < 0
    
    % create and fill tmpHeader struct
    tmpHeader = [];
    
    tmpHeader.fileSizeBytes = fread(FID,1,'uint64');
    tmpHeader.fileType      = char(fread(FID,3,'char')');
    tmpHeader.fileVersion   = fread(FID,1,'char');
    
    % Extrct Header information
    if tmpHeader.fileVersion < 3
        % Event name
        if tmpHeader.fileVersion == 2
            tmpHeader.eventName  = char(fread(FID,4,'char')');
        else
            tmpHeader.eventName  = fliplr(char(fread(FID,4,'char')'));
        end % END IF fileVersion
        
        tmpHeader.channel     = fread(FID, 1, 'uint16'); % Current channel
        tmpHeader.numChan     = fread(FID, 1, 'uint16'); % Number of channels
        tmpHeader.numByteSamp = fread(FID, 1, 'uint16'); % Number of bytes per sample
        reserved              = fread(FID, 1, 'uint16'); % Reserved segment
        
        % Data format in lower four bits
        tmpHeader.dataFormat = allowedFormats{bitand(fread(FID, 1, 'uint8'),7)+1};
        
        % Items used to calculate Fs
        tmpHeader.decimate   = fread(FID, 1, 'uint8');
        tmpHeader.rate       = fread(FID, 1, 'uint16');
        
        % Reserved tags
        reserved = fread(FID, 1, 'uint64');
        reserved = fread(FID, 2, 'uint16');
        
    else
        error(['Unknown file version: ' num2str(tmpHeader.fileVersion)]);
    end % END IF fileVersion
        
    for kk = 1:size(allowedFormats,2)
        
        if strcmp(tmpHeader.dataFormat, allowedFormats{kk})
            numByte = allowedFormatsSize(kk);
        end % END IF
        
    end % END FOR
    % fseek to the start of the segment of interest
%     fseek(FID, floor(time2Extract(1)*Header.Fs) * numByte , 'bof');
%     tmpData = fread(FID, diff(time2Extract)*Header.Fs, ['*' tmpHeader.dataFormat])'; % Read data from file
    tmpData = fread(FID, inf, ['*' tmpHeader.dataFormat])';
    
    varName = ['C', num2str(ii)]; % Name the channel variable
    eval([varName '=tmpData;']); % Set the channel variable to tmpData
    
        % Save Header to file/Create save file
    if headerCount == 0
        Header.info.numSamp = size(tmpData,2);
        Header.info.numChan = length(fileList);
        save(outputPath, 'Header', '-v7.3')
        headerCount = 1;
    end % END headerCount == 0
    
    save(outputPath, varName, '-append') % Save the channel variable
    
    clear(varName);
    
end % END FOR
end % END FUNCTION

% EOF