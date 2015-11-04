% eventIdx - {n}[m] cell with each cell being an event and indicies in the 
%            matrix represent subevents. The first index of the [m] matrix 
%            is the event start.
%
%           Event                         Event
%           ______    ___                  ___    _
%   _______|      |__|   |________________|   |__| |_______
%          Subevent   sub                 sub    sub


function [ eventIdx ] = FindStepEvents(syncData, params )

% Parse Input

% unwrap sync data and transpose.
syncData = syncData(:)';

if ~isfield(params, 'treshold');    threshold = max(syncData)/2; else threshold = params.threshold; end % [Arbitrary units]
if ~isfield(params, 'eventGap');    eventGap  = 4;               else eventGap  = params.eventGap;  end % [Seconds]
if ~isfield(params, 'Fs');          Fs        = 500;             else Fs        = params.Fs;  end       % [Hz]

% normalize to square wave
syncIdx = syncData > threshold;
syncIdx = xor(syncIdx, [0, syncIdx(1:end-1)]); % Find steps
stepIdx = find(syncIdx == 1);

% Remove diff = 1. This is an artifact if the system only records step edges. This happens on XLTec systems.
oneDiffIdx = diff(stepIdx)==1;
stepIdx(oneDiffIdx) = [];

% Remove initial event if it is a partial event
% if diff(stepIdx(1:2)) > eventGap*Fs

% Find event candidates
eventCand = find(diff(stepIdx)>(eventGap*Fs));
numEvents = size(eventCand,2) + 1;

for ii = 1:numEvents
    if ii == numEvents
        eventIdx{ii} = stepIdx(eventCand(ii-1):end);
    else
        if eventCand(ii) == 1
            eventIdx{ii} = stepIdx(1);
        elseif ii == 1
            eventIdx{ii} = stepIdx(1:eventCand(ii));
        else
            eventIdx{ii} = stepIdx(eventCand(ii-1)+1:eventCand(ii));
        end % END IF
    end % END IF
end % END FOR
end % END FUNCTION

% EOF