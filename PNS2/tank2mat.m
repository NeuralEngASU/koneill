% params.eventNames = {'Sync'};
% params.eventChans = 1;

function data = tank2mat( params )

if ~isfield(params, 'eventNames'); eventNames = {'Sync'}; else eventNames = params.eventNames; end
if ~isfield(params, 'eventChans'); eventChans = 1;        else eventChans = params.eventChans; end

% Default Variables
EVENTNAME = {'LFPx','xWav','eNe1'};
CHANNEL   = 1:96;
BLOCKPATH = '';

extList = {'.tev', '.tsq', '.sev'};

% If no path was given, make the user select the block
if isempty(blockPath)
    blockPath = uigetdir('C:\', 'Select Block to Open');
end % END IF

Header = struct('Tank', [], 'TankPath', [], 'Block', [], 'Size', [], 'Type', [],...
    'Code', [], 'Channel', [], 'SortCode', [], 'TimeStart', [],...
    'EventOffset', [], 'DataFormat', [], 'fs', []);

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
TEVFID = fopen(tevFullPath, 'r', 'l');
TSQFID = fopen(tsqFullPath, 'r', 'l');

eventCount = 0;

dataSet = fullfile('C:\CodeRepo\Results\PNS\Data', 'TestFile.h5');

% Loop over eventNames (LFPs, xWAV)
for n = 1%:length(eventName)
    evName = 'Demo';%eventName{2};
    
    evCode = fliplr(double(evName));
    evCode = dec2bin(evCode);
    evCode = [num2str(zeros(4,1)), evCode];
    evCode = reshape(evCode', 1, 32);
    evCode = bin2dec(evCode);
    
    
    
    % Loop over selected channels
    for chan = channel
        
        %h5create(dataSet, ['/C', num2str(chan)], [1, inf])
        
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
            
            %             disp(tsq);
            
            %             if tsq.Code == 1986090872
            %                 disp('Derp')
            %             end % END IF
            
            
            if ~isempty(tsq.Size) && eventCount < 1000
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
                    
                end % END IF
            else
                tsqFlag = false;
            end % END IF
            
            % Measure how many headers are left.
            
            
        end % END WHILE TSQ       
    end % END FOR CHANNELS
end % END FOR EVENT NAME

end % END FUNCTION

% EOF