%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sev2mat.m  (based on SEV2mat.m from TDT)
%   Author: Kevin O'Neill
%   Date: 2015/06/03
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

%% Initialize variables

startTime = ''; % Leave empty if the activeX system works. Else insert the start time of the block you want to extract
timeOfInterst = '10:00:00'; % The time (in real time) of interest. Uses 24hr time

timeBounds = [-5, 5]; % Extract data +/- 5 minutes around the timeOfInterest

sourceDir = ''; % The source directory of the *.sev files (block)
targetDir = ''; % The location where you want to save the extracted data

tank = ''; % The name of the tank
block = ''; % The name of the block with the recording you are interested in

if isempty(sourceDir)
    sourceDir = uigetdir('', 'Select Block Folder');
end % END IF isempty(sourceDir)

if isempty(targetDir)
    targetDir = uigetdir('', 'Select Output Folder');
end % END IF isempty(targetDir)



%% Load tank/block time information

if ~isempty(startTime)
    server = 'Local';
    
    % create TTankX object
    h = figure('Visible', 'off', 'HandleVisibility', 'off');
    TTX = actxcontrol('TTank.X', 'Parent', h);
    
    % connect to server
    if TTX.ConnectServer(SERVER, 'TDT2mat') ~= 1
        close(h)
        error(['Problem connecting to server: ' SERVER])
    end
    
    % open tank
    if TTX.OpenTank(tank, 'R') ~= 1
        TTX.ReleaseServer;
        close(h);
        error(['Problem opening tank: ' tank]);
    end
    
    % select block
    if TTX.SelectBlock(['~' block]) ~= 1
        block_name = TTX.QueryBlockName(0);
        block_ind = 1;
        while strcmp(block_name, '') == 0
            block_ind = block_ind+1;
            block_name = TTX.QueryBlockName(block_ind);
            if strcmp(block_name, block)
                error(['Block found, but problem selecting it: ' block]);
            end
        end
        error(['Block not found: ' block]);
    end
    
    % set info fields
    startTime = TTX.CurBlockStartTime;
    stopTime = TTX.CurBlockStopTime;
    total = stopTime-startTime;
    
    Header.tankPath = TTX.GetTankItem(tank, 'PT');
    Header.date = TTX.FancyTime(startTime, 'Y-O-D');
    Header.startTime = TTX.FancyTime(startTime, 'H:M:S');
    Header.stopTime = TTX.FancyTime(stopTime, 'H:M:S');
    if stopTime > 0
        Header.duration = TTX.FancyTime(total, 'H:M:S');
    else
    Header.duration = -1;
    end
else
    Header.startTime = startTime;
    Header.stopTime = -1;
    Header.duration = -1;
    Header.tankPath = '';
    Header.date = '';
end % END IF

Header.tank = tank;
Header.block = block;

Header.timeBounds = timeBounds;
Header.timeOfInterst = timeOfInterest;

Header.sourceDir = sourceDir;
Header.targetDir = targetDir;

%% Compute number of seconds between startTime and timeOfInterest

timeExpr = '([0-9]+):([0-9]+):([0-9]+)';
startToken = regexp(startTime, timeExpr, 'Tokens');
startSeconds = str2double(startToken{1}{1}) * 3600 + str2double(startToken{1}{2}) * 60 + str2double(startToken{1}{3}) * 1;

toiToken = regexp(timeOfInterest, timeExpr, 'Tokens');
toiSeconds = str2double(toiToken{1}{1}) * 3600 + str2double(toiToken{1}{2}) * 60 + str2double(toiToken{1}{3}) * 1;

timeDiff = toiSeconds - startSeconds;

% If timeOfInterest is smaller than start time (ie: time crosses midnight)
% add 24 hours to the timeOfInterest time.
if timeDiff < 0
    toiSeconds = toiSeconds + 24*3600;
    timeDiff = toiSeconds - startSeconds;
end % END IF

% Compute the time bounds for the segment of interest
time2Extract = [timeDiff + timeBounds(1)*60, timeDiff + timeBounds(2)*60];

% Convert seconds to H:M:S
tmpTime = [toiSeconds + timeBounds(1)*60];
segH = floor(tmpTime/3600);
tmpTime = mod(tmpTime, 3600);
segM = floor(tmpTime/60);
tmpTime = mod(tmpTime, 60);
segS = floor(tmpTime);

segTime = [num2str(segH), ':', num2str(segM), ':', num2str(segS)];

Header.segmentStartTime = segTime;

% Convert seconds H:M:S
tmpTime = [toiSeconds + timeBounds(2)*60];
segH = floor(tmpTime/3600);
tmpTime = mod(tmpTime, 3600);
segM = floor(tmpTime/60);
tmpTime = mod(tmpTime, 60);
segS = floor(tmpTime);

segTime = [num2str(segH), ':', num2str(segM), ':', num2str(segS)];

Header.segmentStopTime = segTime;
%% Extract Header and Test for the length of the sev files

fileList = dir(fullfile(sourceDir, '*.sev'));

if length(fileList) < 1
    warning(['No .sev files found in: ', sourceDir])
end % END IF length(fileList) < 1

% List of allowed formats
allowedFormats = {'single','int32','int16','int8','double','int64'};

% Output the Header once to a .mat file
headerCount = 0;

% Extract the filename for the data.
exprStr = '([A-Za-z0-9\_]+)\_DSP[0-9]\_Ch[0-9]\.[sev]';
fileName = regexp(fileList(1).name, exprStr, 'Tokens');
fileName = fileName{1};

if isempty(targetDir)
    outPath = fullfile(targetDir, [fileName, '.mat']);
end % END IF isempty(targetDir)

filePath = fullfile(sourceDir, fileList(ii).name);
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

numSamp = tmpHeader.fileSizeBytes/ tmpHeader.numByteSamp - 10;
numSec = numSamp / tmpHeader.Fs;

if numSec <= time2Extract(2)
    error('The time-of-interest range is beyond the scope of this file');
end % END IF

Header.Fs = tmpHeader.Fs;
Header.fileVersion = tmpHeader.fileVersion;
Header.numChan = length(fileList);
Header.numSamp = numSamp;
Header.numSampSegment = diff(time2Extract)*Header.Fs;
Header.dataFormat = 'double';
Header.outputPath = tmpHeader.outputPath;

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
    
    % Save Header to file/Create save file
    if headerCount == 0
        save(outPath, 'Header', '-v7.3')
        headerCount = 1;
    end % END headerCount == 0
    
    % fseek to the start of the segment of interest
    fseek(FID, time2Extract(1)*Header.Fs, ['*' tmpHeader.dataFormat]);
    tmpData = fread(FID, diff(time2Extract)*Header.Fs, ['*' tmpHeader.dataFormat])'; % Read data from file
    
    varName = ['C', num2str(ii)]; % Name the channel variable
    eval([varName '=tmpData;']); % Set the channel variable to tmpData
    
    save(outPath, varName, '-append') % Save the channel variable
    
    clear(varName);
    
end % END FOR
% end % END FUNCTION

% EOF