%% How to run this program:
% 
% params.storeNames = {'MicA'  'Sync'  'Mark'};
% params.blockPath = 'D:\Data\HumanECoG\Speech'; % Path to the TDT block.
%   % This item is not needed as the program will uigetdir() you if you
%   % this field is not in the params struct.
% params.targetDir = 'D:\Data\HumanECoG\Speech'; % Path to where the output
%   % .mat file should be located, % This item is not strictly needed as the
%   % program will uigetdir() you if this field is not in the params struct.
% 
% Header = block2mat(params);
%

function Header = block2mat( params )
%% Parse Input
if ~isfield(params, 'storeNames'); storeNames = {'MicA', 'Sync'}; else storeNames = params.storeNames; end 
if ~isfield(params, 'blockPath'); 
    blockPath = uigetdir('C:\', 'Please select block directory'); 
else
    blockPath = params.blockPath;
end

if ~isfield(params, 'targetPath'); 
    targetPath = blockPath; 
else
    targetPath = params.targetPath;
end

% List of extensions
extList = {'.tev', '.tsq', '.sev'};

% Parse pathname
if strcmp(blockPath(end), '\')
    blockPath = blockPath(1:end-1);
end % END IF

pathParts = regexp(blockPath, '\', 'Split');

% Header = struct('tank', [], 'tankPath', [], 'block', [], 'blockPath', [], 'blockSize', [], 'startTimeRaw', []); 
Header = struct('info', [], 'data', []);
Header.info.tank = pathParts{end-1};
Header.info.tankPath = fullfile(pathParts{1:end-1});
Header.info.block = pathParts{end};
Header.info.blockPath = blockPath;
Header.info.blockSize = '';
Header.info.startTime = '';
Header.info.startTimeRaw = [];

for ii = 1:length(storeNames)
   eval(['Header.data.', storeNames{ii}, '.data = [];']);
   eval(['Header.data.', storeNames{ii}, '.Fs = [];']);
end

% Find files
tevList = dir(fullfile(blockPath, ['*', extList{1}]));	% Gets the list tev files
tsqList = dir(fullfile(blockPath, ['*', extList{2}]));  % Gets the list tsq files
% sevList = dir(fullfile(blockPath, ['*', extList{3}]));	% Gets the list sev files

% Look for the names of TEV files
% if ~isempty(tevList)
%     numList = numel(tevList);
%     for i = 1:numList
%         tmpFileList{i} = tevList(i).name;
%     end % END FOR
%     tevList = tmpFileList;
%     clear('tmpFileList');
% end % END IF
% Look for the name of TSQ files
if ~isempty(tevList)
    tevName = tevList(1).name;
end % END IF

% Look for the name of TSQ files
if ~isempty(tsqList)
    tsqName = tsqList(1).name;
end % END IF

% % Look for the name of TSQ files
% if ~isempty(sevList)
%     sevName = sevList(1).name;
% end % END IF

tevFullPath = fullfile(blockPath, tevName);
tsqFullPath = fullfile(blockPath, tsqName);
% sevFullPath = fullfile(blockPath, sevName);

% Open files

eventCount = 0;

% Loop over eventNames (LFPs, xWAV)
for n = 1:length(storeNames)
%     evName = 'Demo';%eventName{2};

    TEVFID = fopen(tevFullPath, 'r', 'l');
    TSQFID = fopen(tsqFullPath, 'r', 'l');

    evName = storeNames{n};
    evCode = fliplr(double(evName));
    evCode = dec2bin(evCode);
    evCode = [num2str(zeros(4,1)), evCode];
    evCode = reshape(evCode', 1, 32);
    evCode = bin2dec(evCode);
    structName = ['Header.data.', evName];
    % Loop over selected channels
    for chan = 1%channel
                
        data = [];
        
        tsqFlag = true;
        eventCount = 0;
        
        % Loop through TSQ
        while tsqFlag
            
            tsq.Size       = fread(TSQFID, [1,1], '*long');
            tsq.Type       = fread(TSQFID, [1,1], '*long');
            tsq.Code       = fread(TSQFID, [1,1], '*long');
            tsq.Channel    = fread(TSQFID, [1,1], '*ushort');
            tsq.SortCode   = fread(TSQFID, [1,1], '*ushort');
            tsq.TimeStart  = fread(TSQFID, [1,1], '*double');
            tsq.EvOffset   = fread(TSQFID, [1,1], '*int64');
            tsq.DataFormat = fread(TSQFID, [1,1], '*long');
            tsq.Fs         = fread(TSQFID, [1,1], '*float');

            if ~isempty(tsq.Size)
                if (tsq.Code == evCode && tsq.Channel == chan)
                    
                    % Move TEV position to the offset specfied in TSQ
                    fseek(TEVFID, tsq.EvOffset, 'bof');
                    
                    % Calculate the number of samples to be gathered
                    numPts = tsq.Size-10;   % Size in longs, minus the header length (long)
                    fetchSize = numPts;     % Number of samples, usually 4 bytes per sample (float)
                    
                    % Read data from file
                    tmpData = fread(TEVFID, [1, fetchSize], '*float');
                    
                    data = [data, tmpData];
                    
                    clear('tmpData')
                    
                    eventCount = eventCount + 1;
                    clc;
                    fprintf('Channel: %d\tEvent: %d\n', chan, eventCount);
                    
                    Fs = tsq.Fs;
                    
                    if eventCount == 1
                        eval([structName, '.Fs = tsq.Fs;']);
                        eval([structName, '.dataFormat = tsq.DataFormat;']);
                        Header.info.startTimeRaw = tsq.TimeStart;
                    end % END IF
                end % END IF
            else
                tsqFlag = false;
            end % END IF
            
        end % END WHILE TSQ        
    end % END FOR CHANNELS
    eval([structName, '.data = data;']);
    
    fclose(TSQFID);
    fclose(TEVFID);
    
end % END FOR EVENT NAME

Header = blocksev2mat(Header, params);

end % END FUNCTION

% EOF