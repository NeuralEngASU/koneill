%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time2Samp
%   Desc: Converts an inputted time to samples (and hours, minutes,
%   seconds)
%   Date: 2015.05.01
%   Author: Kevin O'Neill
%   Principal Invecsigator: Bradley Greger, PhD
%   Lab: Neural Engineering Laboratory, ASU
%
%   [samples] = Time2Samp(timeIn, Fs);
%   [samples, hours, minutes, seconds] = Time2Samp(timeIn, Fs, 'Format', '%H:%M:%S');
%
%   Inputs:
%       timeIn: 'String' of time or a {cell vector} of times. This variable 
%               contains the time(s) that will be converted to samples,
%               hours, minutes, and seconds. Must be in a format like:
%               '1:30:25'. Inputs can be in the form of h:m:s, h:s, m:s
%               
%       Fs:     The sampling rate used to calculate total samples. [Hz]
%
%       Format: (Optional) Entered as a property value: 'Format', 
%               '%H:%M:%S' (default). Formats can be modular, the available
%               formats are:
%
%                   '%H:%M:%S' (default) Long form hours, minutes, seconds.
%                   '%H:%M'    Long form hours and minutes
%                   '%M:%S'    Long form minutes and seconds
%
%                   * Potential Future Functionality * 
%                   '%h:%m:%s' Short form hours, minutes, seconds
%                              Currently same functionality as '%H:%M:%S'.
%
%   Outputs:
%       samples: The number of samples contained in each inputted time 
%               string. Output is a column vector and each entry
%               corrosponds to an inputted time input. [Double]
%
%       hours:  (Optional) The number of hours in each inputted time.
%           Example: '1:30:0' in format '%H:%M:%S'. hours would have a
%                    value of 1.5.
%                    '4572:30' in format '%M:%S'. hours would have a
%                    value of 76.2083.
%
%       minutes: (Optional) The number of minutes in each inputted time.
%               Column vector, each entry corrosponds to each time input.
%
%       seconds: (Optional) The number of seconds in each inputted time.
%               Column vector, each entry corrosponds to each time input.
%               
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ timeSamp, varargout ] = Time2Samp( timeIn, Fs, varargin )

% Default Inputs
FORMAT = '%H:%M:%S';

% Parse Varargin
for ii = 1:2:length(varargin)
    if ~exist(upper(varargin{ii}), 'var')
        fprintf('Unknown option entry: %s\n', varargin{ii})
        return;
    else
        eval([upper(varargin{ii}) '=varargin{ii+1};']);
    end % END IF exist
end % END FOR varargin

formatTime = FORMAT; 

clear FORMAT

%% Count times

% Count the number of times
if iscell(timeIn)
    timeIn = timeIn(:); % Unwrap timeIn
    numTimes = size(timeIn,1); % Count number of times
else
    numTimes = 1;
end % END IF iscell(timeIn)

%% Handle Formats

% Available time codes
timeOptions = {'%H', '%h', '%M', '%m', '%S', '%s'};
timeCodeExpr = '(\%[HMShms])';

% Find time codes in formatTime string
timeCode = regexp(formatTime, timeCodeExpr, 'Match');
numCode = size(timeCode,2); % Number of time codes

% Find which time code options were chosen
idxOptions = zeros(1,6);
for idxCode = 1:numCode
    idxOptions = idxOptions | strcmp(timeOptions, timeCode(idxCode));
end % END FOR numCode

% Define the hms time code as indicies
hmsDef = find(idxOptions);
hmsDef = (hmsDef-1) / 2 + 1;

%===============
% Possible Future functionality
% multiplyOptions = mod(hmsDef,3);
% multiplyOptions = multiplyOptions + mod(multiplyOptions,2) - ~mod(multiplyOptions,2);
% multiplyOptions = multiplyOptions .* (multiplyOptions > 0);
% 
% multiplyOptions = 60.^multiplyOptions;
%================

multiplyOptions = [3600, 60, 1]; % Multiplication factors for time codes.

% Create expression string for each given time code
formatExpr = '';
for idxCode = 1:numCode
    formatExpr = [formatExpr, ':([0-9]+)'];
end % END FOR numCode

% Clean up formatExpr
formatExpr = formatExpr(2:end);

%% Parse Time

% Allocate a matrix of times.
hmsMat = zeros(numTimes,3); % [H,M,S; H,M,S; ...]

% For each inputted time, extract
for numTimeIdx = 1:numTimes
    
    % Extract hms times from string
    tmpMatch = regexp(timeIn{numTimeIdx}, formatExpr, 'Tokens');
    
    % Convert hms times to number
    tmpNum = cellfun(@(x) str2double(x), tmpMatch{1});
    
    % Place converted hms times into their appropriot index 
    hmsMat(numTimeIdx,hmsDef) = tmpNum;
    
end % END FOR numTimes

% Convert times to seconds
timeSec = hmsMat * multiplyOptions'; 

%% Handles Outputs

% Convert to samples
timeSamp = timeSec .* Fs;

% Calculate the number of hours, minutes, seconds
timeTotal = bsxfun(@rdivide, repmat(timeSec,1,3),multiplyOptions);

% Assign outputs
varargout{1} = timeTotal(:,1);
varargout{2} = timeTotal(:,2);
varargout{3} = timeTotal(:,3);

end % END FUNCTION Time2Samp()

% EOF