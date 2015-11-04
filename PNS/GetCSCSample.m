function [ data, timeStamps ] = GetCSCSample( varargin )

timeType = 'sec';

% Variable Initialization
if numel(varargin) < 3
    fprintf('Not enough inputs for GetCSCSample.m\nMinimum Inputs: fileName or '''', channel, timeSpan\n')
elseif numel(varargin) == 3
    fileName = varargin{1};
    channel  = varargin{2};
    timeSpan = varargin{3};
elseif numel(varargin) == 4
    fileName = varargin{1};
    channel  = varargin{2};
    timeSpan = varargin{3};
    timeType = varargin{4};
else
    fprintf('Too many inputs for GetCSCSample.m\nInputs: fileName or '''', timeSpan, timeType\n')
end % END IF

%% Variable Checking

% Time Span Check
if length(timeSpan) ~= 2
    fprintf('timeSpan is incorrectly formatted. Please input data as [beg, end].\n');
elseif (timeSpan(2) - timeSpan(1)) <= 0
    fprintf('timeSpan is incorrectly formatted. Please input data as [beg, end].\n');
end % END IF

% Time type check and format
if strcmp(timeType, 'sec')
    timeSpan = timeSpan * 32000;
elseif strcmp(timeType, 'ts')
    timeSpan = timeSpan;   
else
    fprintf('Input for time type is invalid. Please use ''sec'' or ''ts''\n')
end % END IF

if timeSpan(1) <= 0
    timeSpan(1) = 1;
end;

% Channel number check
channel = ['/CSC', num2str(channel)];

% Data size check
dataLength = timeSpan(2) - timeSpan(1);
dataLength = dataLength * 8;

[userview, ~] = memory;
maxBytes = userview.MaxPossibleArrayBytes;

if dataLength > 0.25*maxBytes
    fprintf('Please use a smaller timeSpan, current size is too large.\nPlease upgrade to 64-bit windows'\n')
end % END IF

%% Data Extraction

data = h5read(fileName, channel, [1, timeSpan(1)], [1, timeSpan(2)-timeSpan(1)]);

% load([fileName(1:end-7), '_CSC.mat'])

timeStamps = [timeSpan(1):timeSpan(2)-1];


end % END FUNCTION

% EOF