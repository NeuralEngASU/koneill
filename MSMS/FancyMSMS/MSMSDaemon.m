%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Title:      MSMSDaemon
%   Desc:       Controls the MSMS object.
%   Author:     Kevin O'Neill (Greger Lab)
%   Date:       August 14, 2012
%
%   Use:        Run the executable MSMSDaemon.exe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% INIT %%%%%

% Params to be passed to m_MSMS
Params.mode           = 'Light';  % Calculation speed (Light/Full)
Params.taskType       = 'Finger'; % Type of task (Finger, Wrist, Arm)
Params.errorLevel     = 'Error';  % Level of Error reported in the Log
Params.scenarioPath   = 'E:\MSMS\MSMS 2.1\Tutorials\FingerTaskDemo'; % Path to the Scenario in MSMS
Params.waitTime       = 0.015;    % Time to wait between UDP packet transmissions to MSMS in seconds. Minimun of 0.01s
Params.loopTime       = 0.5;      % Length of time for the finger to fully open or close in seconds.
Params.Debug.filePath = '';       % Path to output the file created by the EventLog object
Params.Debug.logOut   = false;    % Boolean value to designate whether or not to output the log file.

exitFlag = false;
exprStr = '([A-za-z0-9]+?):([0-1]\.*[0-9]*)';

% Creates the MSMS object and initializes it
m_MSMS = MSMS( );
m_MSMS.Init( Params );


% Clears any UDP object, and creates a new one.
if ~isempty(instrfind('Type','udp')); fclose(instrfind('Type','udp')); end

u = udp('127.0.0.1',4201,'localhost','127.0.0.1','localport',4200);
fopen(u);
if u.BytesAvailable; flushinput(u); flushouput(u); end %clear buffers


%%%%% MAIN %%%%%

% Loops until otherwise told to
while ~exitFlag
    
    inputStr = fscanf(u);
    
    if ~isempty( inputStr )
        
        if strcmp(inputStr, 'Exit')
            exitFlag = true;
        else
            m_MSMS.MoveFingerState( inputStr );
        end % END IF
        
%         inputCell = regexp( inputStr, exprStr, 'Tokens' )';
% 
%         tmpCell = zeros(length(inputCell),2);
%         for k = 1:length(inputCell)
%             tmpCell(k, [1:2]) = str2double(inputCell{k});
%         end % END FOR
%         
%         m_MSMS.MoveFingerState( tmpCell(:,1), tmpCell(:,2) );

    end % END IF
    
    
end % END WHILE

m_MSMS.Destruct( );