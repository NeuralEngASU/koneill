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

function [outPath] = sev2mat(sourceDir, targetDir)

if isempty(sourceDir)
    sourceDir = uigetdir('', 'Select Tank Folder');
end % END IF isempty(sourceDir)

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

    % Save Header to file/Create save file
    if headerCount == 0
        Header = tmpHeader;
        Header = rmfield(Header, 'channel');
        Header.outputPath = outPath;
        
        save(outPath, 'Header', '-v7.3')
        headerCount = 1;        
    end % END headerCount == 0
    
    tmpData = fread(FID, inf, ['*' tmpHeader.dataFormat])'; % Read data from file

    varName = ['C', num2str(ii)]; % Name the channel variable
    eval([varName '=tmpData;']); % Set the channel variable to tmpData
    
    save(outPath, varName, '-append') % Save the channel variable
 
end % END FOR
% end % END FUNCTION

% EOF