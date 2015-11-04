%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Title:      MSMSDaemon
%   Desc:       Controls the MSMS object.
%   Author:     Kevin O'Neill (Greger Lab)
%   Date:       August 14, 2012
%
%   Use:        Run the executable MSMSDaemon.exe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% INIT %%%%%

%%
addpath('E:\MSMS\MSMS 2.1\bin\Matlab\');            % MSMS's MATLAB tools
addpath('F:\CodeRepo\PatientCart\MSMS\Resources\'); % Resources folder location

%%

% Params to be passed to m_MSMS
Params.mode           = 'Light';  % Calculation speed (Light/Full)
Params.taskType       = 3;        % Type of task (ThumbsUp, Fist, FingerExerciser...)
Params.errorLevel     = 'Error';  % Level of Error reported in the Log
Params.scenarioPath   = 'E:\MSMS\MSMS 2.1\Scenarios\MaleRight'; % Path to the Scenario in MSMS
Params.waitTime       = 0.015;    % Time to wait between UDP packet transmissions to MSMS in seconds. Minimun of 0.01s
Params.loopTime       = 0.5;      % Length of time for the finger to fully open or close in seconds.
Params.Debug.filePath = '';       % Path to output the file created by the EventLog object
Params.Debug.logOut   = false;    % Boolean value to designate whether or not to output the log file.

exitFlag = false;
exprStr = '([A-za-z0-9]+?):([0-1]\.*[0-9]*)';


% Creates the MSMS object and initializes it
m_MSMS = MSMS( );
m_MSMS.Init( Params );

% TASK 1 Success
% I - T: 0.3, I: 0.70
% M - T: 0.5, M: 0.75
% R - T: 0.7, R: 0.85
% L - Figure out a better way


%%
% Clears any UDP object, and creates a new one.
if ~isempty(instrfind('Type','udp')); fclose(instrfind('Type','udp')); end

u = udp('127.0.0.1',4201,'localhost','127.0.0.1','localport',4200);
fopen(u);
if u.BytesAvailable; flushinput(u); flushouput(u); end %clear buffers


%%%%% MAIN %%%%%

% Loops until otherwise told to
while ~exitFlag
    
    % Scans for a UDP input
    inputStr = fscanf(u);
    
    if ~isempty( inputStr )
        
        % Detects either a task-switch or an exit string.
        taskSwitch = regexp( inputStr, 'Task=([0-9]+)', 'Tokens' );
        exitDetect = regexp( inputStr, '(Exit)', 'Tokens' );
        
        % If a task-switch was detected, reset MSMS
        % Else if an exit was detected, exit the loop
        % Else, send the string through to MSMS.MoveFingerProp
        if ~isempty( taskSwitch )
            taskSwitch = str2double( taskSwitch{1,1} );
            m_MSMS.TaskSwitch( taskSwitch );
        elseif ~isempty( exitDetect )
            exitFlag = true;
        else
            m_MSMS.MoveFingerProp( inputStr );
        end % END IF
    end % END IF
end % END WHILE

m_MSMS.Destruct( );