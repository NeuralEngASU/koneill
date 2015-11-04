%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MSMSArm < handle % '< handle' means that any instantiation of this class becomes a handle to the object, rather than the obect.
    %MSMSARM Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Properties of MSMSArm
    
    % Only the object itself can set the values to these variables. But
    % 'outside world' can access them.
    properties (SetAccess = private)
        % The 'm_MSMSArm.' before each variable is implicit.
        
        mode = 'Light';
        taskType = 'Finger';
        
        % Contains the structure for the parsed simulationSetup.xml and prostheses.xml
        simInfo   = [];
        
        % Contains the cell for the 'relative Position' matricies for each POI
        relPos    = [];
        
        % Contains a log of erros that happened throughout the object's life
        errorLog = [];
        
        debugOut = true;
        
    end % END PROPERTIES
    
    %% Methods of MSMSArm
    methods
        %% Object Creation/Initialization
        %   Functions to create and initialize the object with the desired
        %   properties
        %
        %   Usage: Line1: obj = CLASSNAME();
        %          Line2: obj.init( params1, params2 );
        
        % Creates the empty object.
        function m_MSMSArm = MSMSArm()
            % Leave Empty
        end % END FUNCTION
        
        % Initializes the object
        function init( m_MSMSArm, params )
            
            if nargin < 2
                errorMsg = dbstack( );
                m_MSMSArm.Debug( errorMsg, '_NoParamsInput_');
            end % END IF
            
            % Defaults to a faster calculation if no other is selected.
            if isfield(params, 'mode')
                m_MSMSArm.mode = 'Light';
                
                % Not setting params.mode is considered an error and will be
                % recorded accordingly.
                errorMsg = dbstack( );
                m_MSMSArm.Debug( errorMsg, '_NoModeSelected_');
            else
                m_MSMSArm.mode = char(params.mode);
            end% END IF
            
            % Defaults to Finger task if no other is selected.
            if isfield(params, 'taskType')
                m_MSMSArm.taskType = 'Finger';
                
                % Not sending params.taskType is considered an error and will be
                % recorded accordingly
                errorMsg = dbstack( );
                m_MSMSArm.Debug( errorMsg, '_NoTaskTypeSelected_');
            else
                m_MSMSArm.taskType = char(params.taskType);
            end% END IF
            
            if isfield(params, 'controlType')
                m_MSMSArm.taskType = 'State';
                
                % Not sending params.controlType is considered an error and will be
                % recorded accordingly
                errorMsg = dbstack( );
                m_MSMSArm.Debug( errorMsg, '_NoControlTypeSelected_');
            else
                m_MSMSArm.controlType = params.controlType;
            end
            
            % Parses the simulationSetup.xml and prostheses.xml files from
            % the MSMS Scenario
            m_MSMSArm.SimInfo = ParseXML('E:\MSMS\MSMS 2.1\Tutorials\FingerTask');
            
            % Parses and organizes data to form the relPos matrix.
            m_MSMSArm.ParseRelPos( );
            
            % Sets up the scene with the specified taskType
            m_MSMSArm.SetupScene( );
            
        end % END INIT
        
        % END Object Creation/Initialization
        %% Setters
        %   Functions that allow the 'Outside World' to set the variables
        %   of m_MSMSarm manually.
        %
        %   Usage: obj.SetThisVariable( NewValue )
        
        % END Setters
        %% Getters
        %   Function that allow the 'Outside World' to get the value of the
        %   variables within m_MSMSArm.
        %
        %   Usage: value = obj.getThisVariable()
        
        % Gets the current list of errors
        function [ errorLog ] = getErrorList( m_MSMSArm )
            errorLog = m_MSMSArm.errorLog;
        end % END FUNCTION
        
        % Gets the current 'relative Position' cell (relPos)
        function [ relPos ] = getRelPos( m_MSMSArm )
            relPos = m_MSMSArm.relPos;
        end % END FUNCTION
        
        
        % END Getters
        %% Working Bits
        %   Functions that do all of the grunt work. (Calculation,
        %   analysis, data aggregation, and everything else more
        %   complicated that get or set.
        %
        %   Usage:  obj.WorkingBit1( inputs );
        %           val = obj.WorkingBit2( inputs );
        
        % Sets up the current MSMS scene
        function SetupScene( m_MSMSArm, taskType )
            
            if nargin == 2
               m_MSMSArm.taskType = taskType;
            end % END IF
            
            switch m_MSMSArm.taskType
                case 'Finger'
                    % Perhaps a 2 column cell, col1 for feature name, col2
                    % for feature count.
                    m_MSMSArm.featureList = {...
                        'ShoulderRoll'; 'ElbowPitch'; 'RotateWrist';...
                        'PlaceCamera' ;...
                        'PlaceSphere1'; 'PlaceSphere2'; 'PlaceShpere3';...
                        'ColorSphere1'; 'ColorSphere2'; 'ColorSphere3';...
                        'HideSphere1' ; 'HideSphere2' ; 'HideSphere3' ;...
                        'ZeroFingers'};
                otherwise
                    errorMsg = dbstack( );
                    m_MSMSArm.Debug( errorMsg, ['_NoCaseForTaskTypeInput_', m_MSMSArm.taskType]);
                    
            end % END SWITCH
            
            % Creates a UDP packet with the selected featureList
            m_MSMSArm.CreateUDP;
            
%%%%%%%%%%%%% Will need to think of a way to parse the featureList

            m_MSMSArm.ShoulderRoll;
            m_MSMSArm.ElbowRoll;
            m_MSMSArm.RotateWrist;
            m_MSMSArm.PlaceCamera;
            m_MSMSArm.PlaceSphere;
            m_MSMSArm.ColorSphere;
            m_MSMSArm.HideSphere;
            m_MSMSArm.ZeroFingers;
            
%%%%%%%%%%%%%

            m_MSMSArm.SendUDP;
            
        end % END FUNCTION
        
        % Handles the retreval of any errors and stores them in an
        % ErrorList
        function Debug( m_MSMSArm, tmpErrorMsg, tmpErrorStr)
            
            % Stores the current error and reason into errorList
            m_MSMSArm.errorLog{end,:} = [ tmpErrorMsg.file, tmpErroMsg.name, tmpErrorMsg.line, tmpErrorStr ];
                  
            % Used to decide what to do with the errorList
            intDescision = m_MSMSArm.debugOut + m_MSMSArm.debugFileDNE;
            
            switch intDescision
                case 1
                    % Get time (StampTime?), Make file, output to it
                case 3
                    % Output to file
                otherwise
                    % Return
            end % END SWITCH
        end % END FUNCTION
        % END WORKING BITS
    end % END METHODS
    
end % END CLASS

