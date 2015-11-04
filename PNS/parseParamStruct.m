function  trialStruct = parseParamStruct( paramStruct )
% Parses the paramStruct structure from parseNEV_FingerPress_Kevin into a
% better form that is comparible with the koneill library.
%
% Returns a structure that has lenght(trial types), each 
%
% Kevin O'Neill
% 20130514

% Initializes output structure
trialStruct = [];

% Does not loop, intended to that the code is collapsable and look
% seperate.
% This section finds the field names to be used in trialStruct
for k = 1:1
    
    %%% Extract marker names from paramStruct
    %Get field names
    sNames = fieldnames(paramStruct(k));
    
    % Find fields which have 'Frame' in the name
    tmpNames = regexp(sNames(1:[end-100]), '(Frame)');
    
    % Create a boolean array to store indicies of interest.
    boolNames = zeros(length(tmpNames), 1);
    
    % loop over indicies in tmpNames and find the empty matricies
    for j = 1:length(tmpNames)
        boolNames(j) = isempty(tmpNames{j});
    end % END FOR
    
    % Convert to a logical column vector
    boolNames = logical(boolNames');
    
    % Pull field names from sNames
    markerNames = sNames(boolNames);
    
    % Seperate marker names
    trialTypes = [regexp(markerNames, '([A-Za-z]+)([0-9]+)', 'Tokens')];
    for j = 1:length(trialTypes)
        trialTypes(j, 1:2) = trialTypes{j}{1};
    end % END FOR
    
    % Sorts trial types by trial number and makes the trial number a valid
    % variable name.
    trialTypes = sortrows(trialTypes,2);
    
    % Find unique trial names
    uniqueTrial = cellfun(@(x) [x], unique(trialTypes(:,2)), 'uniformoutput',false);
    uniqueTrialLen = length(uniqueTrial);
    
    % Create trialStruct structure, one matrix for each trial. [trial#,
    % acquiredTS, thresholdTS, onTS, timeoutTS]
    for j = 1:uniqueTrialLen
        eval(['trialStruct.',strcat('T', uniqueTrial{j}),' = [];']);
    end % END IF
    
end % END FOR

% loops over all trial types
for k = 1:uniqueTrialLen
    % resets currentTrial value
    currentTrial = [];
    % Loops over all trials
    for j = 1:length(paramStruct)
        % checks if the trial is of a trial type AND that the trial is not
        % a 'reset'
        if ~isempty(eval(['paramStruct(j).On', uniqueTrial{k}, 'TS'])) && ~strcmp(paramStruct(j).SS_TrialType, '''reset''')
            % initializes the trial's values
            currentTrial(end+1, 1:5) = zeros(1,5);
            
            % stores the trial value
            currentTrial(end, 1) = str2num(paramStruct(j).SS_TrialCount);
            
            % Stores the Acquired time stamp
            if ~isempty(eval(['paramStruct(j).Acquired', uniqueTrial{k}, 'TS']))
                currentTrial(end, 2) = eval(['paramStruct(j).Acquired', uniqueTrial{k}, 'TS']);
            else
                currentTrial(end, 2) = [NaN];
            end% END IF
            
            % Stores the CrosedThreshold time stamp
            if ~isempty(eval(['paramStruct(j).CrossedThreshold', uniqueTrial{k}, 'TS']))
                currentTrial(end, 3) = eval(['paramStruct(j).CrossedThreshold', uniqueTrial{k}, 'TS']);
            else
                currentTrial(end, 3) = [NaN];
            end% END IF
            
            % Stores the On time stamp
            if ~isempty(eval(['paramStruct(j).On', uniqueTrial{k}, 'TS']))
                currentTrial(end, 4) = eval(['paramStruct(j).On', uniqueTrial{k}, 'TS']);
            else
                currentTrial(end, 4) = [NaN];
            end% END IF
            
            % Checks if there was a timeout field BITWISE AND 
            % if there is a field, if it has a value
            % then stores the Timeout time stamp
            if isfield(paramStruct(j), ['Timeout', uniqueTrial{k}, 'TS']) & ~isempty(eval(['paramStruct(j).Timeout', uniqueTrial{k}, 'TS']))
                currentTrial(end, 5) = eval(['paramStruct(j).Timeout', uniqueTrial{k}, 'TS']);
            else
                currentTrial(end, 5) = [NaN];
            end% END IF
        end % END IF
    end % END FOR
    
    % Stores the trial information in trialStruct
    eval(['trialStruct.',strcat('T', uniqueTrial{k}),' = currentTrial;']);
end % END FOR
end % END FUNCTION

% EOF