%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	FindLaserEvents
%		Kevin O'Neill
%		Mario Capecchi Lab
%		20131118
%		v0.1
%		PI: Naveen Nagarajan
%		
%	Inputs:
%       CSCData: Data loaded from the expYYYY-MM-DD_hh-mm-ss_CSC.mat file
%											(nsc2mat output)
%
%       threshold: The aribrary value by which you want to select for laser
%                  data. The recommended and default value is:
%                           1/2 * max(laserData)
%
%	Outputs:
%       laserEvents: A two column matrix in which the first column contains
%                    the laser on points and the second column contains the
%                    paired laser off indicies.
%
%	To Use:
%       Run the function with the correct inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ trialEvents, laserEvents ] = FindLaserEvents( CSCData, varargin )

% Loads timestamps
laserTime = CSCData.timeStamps;

% Loads laserData from HDF5 file (channel 17).
tempChanName = ['CSC',num2str(17)];
laserData = h5read(CSCData.hdf5FileName, ['/', tempChanName],...
    1,length(laserTime));

if ~isempty(varargin)
    threshold = varargin{1};
else
    threshold = 0.5 * max(laserData);
end % END IF

% Find events
laserData = laserData(:)'; % Unwrap data
laserIdx = laserData > threshold; % Convert to square wave
laserIdx = xor(laserIdx, [0, laserIdx(1:end-1)]); % Find transitions

% laserOffIdx = laserOffIdx < threshold;
% laserOffIdx = xor(laserOffIdx, [laserOffIdx(2:end),0];

% Find TimeStamps
laserEvents = laserTime(laserIdx);

% Reshape to [nx2] matrix
laserEvents = reshape(laserEvents,2, length(laserEvents)/2)';

trialIdx = find(diff(laserEvents(:,1))>1e6); % Find gaps of around 1 second

trialIdxStart = [1, trialIdx+1]; % Shift index to start of trial. 
trialIdxEnd = [trialIdx, size(laserEvents, 1)]; % Shift index to end of trial. 

% Creates the start and end pairs for trials.
trialEvents = [laserEvents(trialIdxStart, 1), laserEvents(trialIdxEnd, 2)];

laserPower = round((laserEvents(:,2)-laserEvents(:,1))/2 + laserEvents(:,1));
laserPower = laserData(find(laserTime == laserPower));

end % END FUNCTION

% EOF