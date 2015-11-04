%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Title:      MSMS
%   Desc:       Controls MSMS via UDP packet building
%   Author:     Kevin O'Neill (Greger Lab)
%   Date:       August 14, 2012
%
%   Use:        m_MSMSobject = MSMS();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% '< handle' means that any instantiation of this class becomes a handle to
% the object, rather than the obect.
classdef MSMS < handle 
    %MSMS Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Properties of MSMS
    
    % Only the object itself can set and get these variables.
    properties
        % The 'm_MSMS.' before each variable is implicit.
        
        % Default params
        mode       = '0';
        taskType   = '0';
        errorLevel = '0';
        
        % Contains the structure for the parsed simulationSetup.xml and prostheses.xml
        SimInfo   = [];
        
        % Contains the cell for the 'relative Position' matricies for each POI
        relPos    = [];
        
        % Path to the scenario for MSMS
        scenarioPath = '';

        % The UDP packet that will communicate with MSMS
        udpPacket = [];
        
        % List of features
        featureList = [];
        
        % Success flags for each POI
        successFlag = false(6);
        
        % States for each POI
        stateFlag = zeros(6);
        
        % Expression for Regular Expressions used when inturpetting the 
        % input string.
        exprStr = '([A-za-z0-9]+?)=([0-1]\.*[0-9]*)';
        
        %%%% Only used during MoveFingerState %%%%
        
        % Time to wait between loops
        waitTime = 0.015; % Miliseconds
        
        % Time the state movement loop takes to complete
        loopTime = 0.5;
        
        % State movement steps 
        stateSteps = 1;
        
        % Matrices to hold maximum/minimum joint angles for DIP/PIP/Knuckle/+ of each finger
        jointAngleMax = zeros(5,4);
        jointAngleMin = zeros(5,4);
        
        % Matrix to hold the joint angles for each joint as calculated by
        % AngleFun
        jointAngle = zeros(5,3);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Colors avaliable to spheres enough POIs 1-6 are also colors 1-6
        colors = [ 1, 0, 0;... % Red
                   0, 1, 0;... % Green
                   0, 0, 1;... % Blue
                   0, 1, 1;... % Cyan
                   1, 0, 1;... % Magenta
                   1, 1, 0;... % Yellow
                   1, 1, 1;... % White
                   0, 0, 0];   % Black
            
        % Sphere Positions
%         spherePos = zeros(6,3);
        % Left
%         spherePos = [ 0.5400, -0.1400, -0.0260;...
%                       0.5750, -0.1435, -0.0110;...
%                       0.5850, -0.1610,     0.0;...
%                          0.8,     0.0,     0.0;...
%                          0.0,     0.0,     0.0;...
%                          0.0,     0.0,     0.0];  
        % Right
        spherePos = [ 0.5800, -0.1575,  0.0120;...
                      0.5850, -0.1530, -0.0175;...
                      0.5900, -0.1725, -0.0300;...
                      0.5830, -0.1910, -0.0400;...
                      0.5600, -0.2055, -0.0520;...
                         0.0,     0.0,     0.0];  
        
        % Position and rotation for the camera
%         cameraPosRot = [ 0.50, (-30 / 360) * (2*pi) + (2*pi);...
%                          0.09, 0.0;...
%                          0.45, 0.0];

        cameraPosRot = [ 0.50, (150 / 360) * (2*pi);...
                        0.075,                  0.0;...
                        -0.45, (180 / 360) * (2*pi)];
        %%%%%%%%%%% Keeping indent in check
        
                % Male and Female models
        Male   = 0;
        Female = 1;
        
        % Left and Right arms
        Left   = 0;
        Right  = 1;
        
        % Quick and long calculation settings
        Light  = 0;
        Full   = 1;
        
        % Debug Log Levels
        Error   = 0;
        Warning = 1;
        Info    = 2;
        
        % Points of Interest
        Thumb  = 1;
        Index  = 2;
        Middle = 3;
        Ring   = 4;
        Little = 5;
        Palm   = 6;
        
        % Spheres
        Sphere1 = 1;
        Sphere2 = 2;
        Sphere3 = 3;
        Sphere4 = 4;
        Sphere5 = 5;
        Sphere6 = 6;
        
        % Joint Information
        POIOffset     = [ 1, 1];
        DIP           = [ 2, 2];
        PIP           = [ 3, 3];
        Knuckle       = [ 4, 4];
        Side          = [ 5, 5];
        WristYaw      = [ 6, 0];
        WristPitch    = [ 7, 0];
        WristRoll     = [ 8, 0];
        ElbowPitch    = [ 9, 0];
        ShoulderRoll2 = [10, 0];
        ShoulderPitch = [11, 0];
        ShoulderRoll  = [12, 0];
        Shoulder1     = [13, 0];
        GlobalXRot    = [14, 6];
        GlobalYRot    = [15, 7];
        GlobalZRot    = [16, 8];
        GlobalOffset  = [17, 9]; 
        
        % Tasks
%         FingerThumb     = 1;
%         FingerPalm      = 2;
%         FingerExerciser = 3;
%         Fist            = 4;
%         ThumbsUp        = 5;
        
    end % END PROPERTIES
    
    %% Methods of MSMS
    methods
        %% Constructor [Object Creation/Initialization]
        %   Functions to create and initialize the object with the desired
        %   properties
        %
        %   Usage: Line1: obj = CLASSNAME();
        %          Line2: obj.init( params1, params2 );
        
        % Creates the empty object.
        function m_MSMS = MSMS()
            % Leave Empty
        end % END FUNCTION
        
        % Initializes the object
        function Init( m_MSMS, Params )    
            
%             m_MSMS.spherePos(:,:) = [ 0.5400, -0.1400, -0.0260;...
%                                       0.5750, -0.1435, -0.0110;...
%                                       0.5850, -0.1610,     0.0;...
%                                          0.8,     0.0,     0.0;...
%                                          0.0,     0.0,     0.0;...
%                                          0.0,     0.0,     0.0];        
            
            % Parses Params
            m_MSMS.mode       = eval(['m_MSMS.', Params.mode, ';']);
            m_MSMS.taskType   = Params.taskType;%eval(['m_MSMS.', Params.taskType, ';']);
            m_MSMS.errorLevel = eval(['m_MSMS.', Params.errorLevel, ';']);
            m_MSMS.waitTime   = Params.waitTime;
            m_MSMS.loopTime   = Params.loopTime;
            m_MSMS.stateSteps = m_MSMS.loopTime/m_MSMS.waitTime;
            
            % Parses the simulationSetup.xml and prostheses.xml files from
            % the MSMS Scenario
            m_MSMS.SimInfo = ParseXML( Params.scenarioPath );
            
            % Sets up the scene with the specified taskType
            m_MSMS.SetupScene( );
            
        end % END INIT
        
        % END Constructor [Object Creation/Initialization]
        %% Destructor
        
        % Destructs the m_MSMS object
        function Destruct( m_MSMS )
            
            % Destruct/output/save anything needed here
            
        end % END FUNCTION
        
        % END Destructor
        %% Setters
        %   Functions that allow the 'Outside World' to set the variables
        %   of m_MSMS manually.
        %
        %   Usage: obj.SetThisVariable( NewValue )
        
        % END Setters
        %% Getters
        %   Function that allow the 'Outside World' to get the value of the
        %   variables within m_MSMS.
        %
        %   Usage: value = obj.getThisVariable()
                
        % END Getters
        %% Working Bits
        %   Functions that do all of the grunt work. (Calculation,
        %   analysis, data aggregation, and everything else more
        %   complicated that get or set.
        %
        %   Usage:  obj.WorkingBit1( inputs );
        %           val = obj.WorkingBit2( inputs );
        
        % Sets up the current MSMS scene
        function SetupScene( m_MSMS )
            
            m_MSMS.jointAngleMin = zeros(5,4);
            
            switch m_MSMS.taskType
                case 1 % Finger-to-Thumb
                    m_MSMS.featureList = {'ShoulderRoll', 1;   
                                            'ElbowPitch', 1;
                                           'RotateWrist', 1;
                                           'PlaceCamera', 1;
                                           'PlaceSphere', 6;
                                           'ColorSphere', 6;
                                            'HideSphere', 6;
                                           'ZeroFingers', 16};
                    
                    % Maximum joint angles figners can reach
                    % Joint columns:
                    %   DIP, PIP, Knuckle, Knuckle-Side
                    %
                    % Joint Rows:
                    %   Thumb
                    %   Index
                    %   Middle
                    %   ...
                    m_MSMS.jointAngleMax = deg2rad( [ 40, 30, -40, 68;...
                                                       0, 90,  90,  0;...
                                                       0, 90,  90,  0;...
                                                       0, 90,  90,  0;...
                                                       0, 90,  90,  0]);
                                                   
                    % Sets current joint angle !!!!DATA!!!! to 0
                    m_MSMS.AngleFun( [0,0,0,0,0] );
                    
                    % Creates a UDP packet with the selected featureList
                    m_MSMS.CreateUDP( sum([m_MSMS.featureList{1:end,2}]) );

                    %%%%%%%%%%%%%

                    % Sets joint angles and opject properties.
                    m_MSMS.MoveJoint( {'ShoulderRoll', 'ElbowPitch', 'WristRoll' },...
                       [(45/360)*(2*pi), (45/360)*(2*pi), (-30/360)*(2*pi) + 2*pi],...
                       [1], false)

                    m_MSMS.PlaceCamera( false );
                    m_MSMS.PlaceSphere( [1:6], false );
                    m_MSMS.ColorSphere( [1:6], 127 , false);
                    m_MSMS.HideSphere(  [1:6], ones(1,6).*100, false);

                    m_MSMS.MoveJoint( { 'ThumbDIP', 'ThumbKnuckle',     'ThumbSide', 'ThumbMeta'}, [0,0, m_MSMS.jointAngleMax(1,3), 0], [1], false );
                    m_MSMS.MoveJoint( { 'IndexDIP',     'IndexPIP',  'IndexKnuckle'}, m_MSMS.jointAngle(2,1:3), [2], false );
                    m_MSMS.MoveJoint( {'MiddleDIP',    'MiddlePIP', 'MiddleKnuckle'}, m_MSMS.jointAngle(3,1:3), [3], false );
                    m_MSMS.MoveJoint( {  'RingDIP',      'RingPIP',   'RingKnuckle'}, m_MSMS.jointAngle(4,1:3), [4], false );
                    m_MSMS.MoveJoint( {'LittleDIP',    'LittlePIP', 'LittleKnuckle'}, m_MSMS.jointAngle(5,1:3), [5], false );

                    %%%%%%%%%%%%%

                    m_MSMS.SendUDP;
                    
                case 2 % Finger-to-Palm
                    m_MSMS.featureList = {'ShoulderRoll', 1;   
                                            'ElbowPitch', 1;
                                           'RotateWrist', 1;
                                           'PlaceCamera', 1;
                                           'PlaceSphere', 6;
                                           'ColorSphere', 6;
                                            'HideSphere', 6;
                                           'ZeroFingers', 16};
                                       
                    m_MSMS.jointAngleMax = deg2rad( [ 40, 30, -40, 68;...
                                                       0, 90,  90,  0;...
                                                       0, 90,  90,  0;...
                                                       0, 90,  90,  0;...
                                                       0, 90,  90,  0]);
                    
                    % Sets current joint angle !!!!DATA!!!! to 0
                    m_MSMS.AngleFun( [0,0,0,0,0] );
                    
                    % Creates a UDP packet with the selected featureList
                    m_MSMS.CreateUDP( sum([m_MSMS.featureList{1:end,2}]) );

                    %%%%%%%%%%%%%
                    
                    % Sets joint angles and opject properties.
                    m_MSMS.MoveJoint( {'ShoulderRoll', 'ElbowPitch', 'WristRoll' },...
                       [(45/360)*(2*pi), (45/360)*(2*pi), (-30/360)*(2*pi) + 2*pi],...
                       [1], false)

                    m_MSMS.PlaceCamera( false );
                    m_MSMS.PlaceSphere( [1:6], false );
                    m_MSMS.ColorSphere( [1:6], 127 , false);
                    m_MSMS.HideSphere(  [1:6], ones(1,6).*100, false);

                    m_MSMS.MoveJoint( { 'ThumbDIP', 'ThumbKnuckle',     'ThumbSide', 'ThumbMeta'}, [0,0, m_MSMS.jointAngleMax(1,3), 0], [1], false );
                    m_MSMS.MoveJoint( { 'IndexDIP',     'IndexPIP',  'IndexKnuckle'}, m_MSMS.jointAngle(2,1:3), [2], false );
                    m_MSMS.MoveJoint( {'MiddleDIP',    'MiddlePIP', 'MiddleKnuckle'}, m_MSMS.jointAngle(3,1:3), [3], false );
                    m_MSMS.MoveJoint( {  'RingDIP',      'RingPIP',   'RingKnuckle'}, m_MSMS.jointAngle(4,1:3), [4], false );
                    m_MSMS.MoveJoint( {'LittleDIP',    'LittlePIP', 'LittleKnuckle'}, m_MSMS.jointAngle(5,1:3), [5], false );

                    %%%%%%%%%%%%%

                    m_MSMS.SendUDP;
                    
                case 3 % Finger Exerciser

                    m_MSMS.featureList = {'ShoulderRoll', 1;   
                                            'ElbowPitch', 1;
                                           'RotateWrist', 1;
                                           'PlaceCamera', 1;
                                           'PlaceSphere', 6;
                                           'ColorSphere', 6;
                                            'HideSphere', 6;
                                           'ZeroFingers', 15};
                                       
                    m_MSMS.jointAngleMin = deg2rad( [ 45, 20,  0, 0;...
                                                      10, 70, 30, 0;...
                                                      10, 70, 30, 0;...
                                                      10, 70, 30, 0;...
                                                      10, 70, 30, 0]);
                                                  
                                                  
                    m_MSMS.jointAngleMax = deg2rad( [ 90, 45,  0, 0;...
                                                      45, 90, 40, 0;...
                                                      45, 90, 40, 0;...
                                                      45, 90, 40, 0;...
                                                      45, 90, 40, 0]);
                                                  
                    m_MSMS.AngleFun( [0,0,0,0,0] );
                    
                    % Creates a UDP packet with the selected featureList
                    m_MSMS.CreateUDP( sum([m_MSMS.featureList{1:end,2}]) );

                    %%%%%%%%%%%%%

                    m_MSMS.MoveJoint( {'ShoulderRoll', 'ElbowPitch', 'WristRoll' },...
                       [(45/360)*(2*pi), (45/360)*(2*pi), (-30/360)*(2*pi) + 2*pi],...
                       [1], false)

                    m_MSMS.PlaceCamera( false );
                    m_MSMS.PlaceSphere( [1:6], false );
                    m_MSMS.ColorSphere( [1:6], 127 , false);
                    m_MSMS.HideSphere(  [1:6], ones(1,6).*100, false);

                    m_MSMS.MoveJoint( { 'ThumbDIP', 'ThumbKnuckle',     'ThumbSide', 'ThumbMeta'}, m_MSMS.jointAngle(1,1:4), [1], false );
                    m_MSMS.MoveJoint( { 'IndexDIP',     'IndexPIP',  'IndexKnuckle'}, m_MSMS.jointAngle(2,1:3), [2], false );
                    m_MSMS.MoveJoint( {'MiddleDIP',    'MiddlePIP', 'MiddleKnuckle'}, m_MSMS.jointAngle(3,1:3), [3], false );
                    m_MSMS.MoveJoint( {  'RingDIP',      'RingPIP',   'RingKnuckle'}, m_MSMS.jointAngle(4,1:3), [4], false );
                    m_MSMS.MoveJoint( {'LittleDIP',    'LittlePIP', 'LittleKnuckle'}, m_MSMS.jointAngle(5,1:3), [5], false );

                    %%%%%%%%%%%%%

                    m_MSMS.SendUDP;
                
                case 4 % Make-A-Fist
                    m_MSMS.featureList = {'ShoulderRoll', 1;   
                                            'ElbowPitch', 1;
                                           'RotateWrist', 1;
                                           'PlaceCamera', 1;
                                           'PlaceSphere', 6;
                                           'ColorSphere', 6;
                                            'HideSphere', 6;
                                           'ZeroFingers', 16};
                                       
                    m_MSMS.jointAngleMax = deg2rad( [  0, 20,  0, 0;...
                                                      90, 90, 90, 0;...
                                                      90, 90, 90, 0;...
                                                      90, 90, 90, 0;...
                                                      90, 90, 90, 0]);
                                                  
                    m_MSMS.AngleFun( [0,0,0,0,0] );
                    
                    % Creates a UDP packet with the selected featureList
                    m_MSMS.CreateUDP( sum([m_MSMS.featureList{1:end,2}]) );

                    %%%%%%%%%%%%%

                    m_MSMS.MoveJoint( {'ShoulderRoll', 'ElbowPitch', 'WristRoll' },...
                       [(45/360)*(2*pi), (45/360)*(2*pi), (-30/360)*(2*pi) + 2*pi],...
                       [1], false)

                    m_MSMS.PlaceCamera( false );
                    m_MSMS.PlaceSphere( [1:6], false );
                    m_MSMS.ColorSphere( [1:6], 127 , false);
                    m_MSMS.HideSphere(  [1:6], ones(1,6).*100, false);

                    m_MSMS.MoveJoint( { 'ThumbDIP', 'ThumbKnuckle',     'ThumbSide', 'ThumbMeta'}, m_MSMS.jointAngle(1,1:4), [1], false );
                    m_MSMS.MoveJoint( { 'IndexDIP',     'IndexPIP',  'IndexKnuckle'}, m_MSMS.jointAngle(2,1:3), [2], false );
                    m_MSMS.MoveJoint( {'MiddleDIP',    'MiddlePIP', 'MiddleKnuckle'}, m_MSMS.jointAngle(3,1:3), [3], false );
                    m_MSMS.MoveJoint( {  'RingDIP',      'RingPIP',   'RingKnuckle'}, m_MSMS.jointAngle(4,1:3), [4], false );
                    m_MSMS.MoveJoint( {'LittleDIP',    'LittlePIP', 'LittleKnuckle'}, m_MSMS.jointAngle(5,1:3), [5], false );

                    %%%%%%%%%%%%%

                    m_MSMS.SendUDP;
                    
                case 5 % Thumbs Up
                    m_MSMS.featureList = {'ShoulderRoll', 1;   
                                            'ElbowPitch', 1;
                                           'RotateWrist', 1;
                                           'PlaceCamera', 1;
                                           'PlaceSphere', 6;
                                           'ColorSphere', 6;
                                            'HideSphere', 6;
                                           'ZeroFingers', 16};
                                       
                    m_MSMS.jointAngleMin = deg2rad( [  0, -80,  0,  0;...
                                                      90,  90, 90, 0;...
                                                      90,  90, 90, 0;...
                                                      90,  90, 90, 0;...
                                                      90,  90, 90, 0]); 
                                       
                    m_MSMS.jointAngleMax = deg2rad( [  0, 20,  0, 0;...
                                                      90, 90, 90, 0;...
                                                      90, 90, 90, 0;...
                                                      90, 90, 90, 0;...
                                                      90, 90, 90, 0]);
                                                  
                    m_MSMS.AngleFun( [0,0,0,0,0] );
                    
                    % Creates a UDP packet with the selected featureList
                    m_MSMS.CreateUDP( sum([m_MSMS.featureList{1:end,2}]) );

                    %%%%%%%%%%%%%

                    m_MSMS.MoveJoint( {'ShoulderRoll', 'ElbowPitch', 'WristRoll' },...
                       [(45/360)*(2*pi), (45/360)*(2*pi), (-30/360)*(2*pi) + 2*pi],...
                       [1], false)

                    m_MSMS.PlaceCamera( false );
                    m_MSMS.PlaceSphere( [1:6], false );
                    m_MSMS.ColorSphere( [1:6], 127 , false);
                    m_MSMS.HideSphere(  [1:6], ones(1,6).*100, false);

                    m_MSMS.MoveJoint( { 'ThumbDIP', 'ThumbKnuckle',     'ThumbSide'}, m_MSMS.jointAngle(1,1:3), [1], false );
                    m_MSMS.MoveJoint( { 'IndexDIP',     'IndexPIP',  'IndexKnuckle'}, m_MSMS.jointAngle(2,1:3), [2], false );
                    m_MSMS.MoveJoint( {'MiddleDIP',    'MiddlePIP', 'MiddleKnuckle'}, m_MSMS.jointAngle(3,1:3), [3], false );
                    m_MSMS.MoveJoint( {  'RingDIP',      'RingPIP',   'RingKnuckle'}, m_MSMS.jointAngle(4,1:3), [4], false );
                    m_MSMS.MoveJoint( {'LittleDIP',    'LittlePIP', 'LittleKnuckle'}, m_MSMS.jointAngle(5,1:3), [5], false );

                    %%%%%%%%%%%%%

                    m_MSMS.SendUDP;
                                                                                  
                otherwise                    
            end % END SWITCH            
        end % END FUNCTION
        
        % Creates the UDP Packet
        function CreateUDP( m_MSMS, featureCount )
            m_MSMS.udpPacket = []; % Empty packet
            m_MSMS.udpPacket = [ m_MSMS.udpPacket; int8(featureCount) ];
        end % END FUNCTION
        
        % Sends the UDP packet to MSMS
        function SendUDP( m_MSMS )
            judp( 'send', 11114, '127.0.0.1', typecast( m_MSMS.udpPacket, 'int8' ) );
        end % END FUNCTION
        
        % Places the given Spheres
        function PlaceSphere( m_MSMS, targets, singleton )
            
            % If no input for singleton is given, assume true
            if nargin < 3
                singleton = true;
            end % END IF
            
            % If function was called as a singleton, create and send the
            % UDP packet
            if singleton
                m_MSMS.CreateUDP( length(targets) * 1 );
                % DEBUG
            else
                % DEBUG out that func in not a singleton
            end % END IF
            
            % Iterates over the number of targets given.
            for k=1:length(targets)
                
                % Parses Sphere info
                eval(['compType   = m_MSMS.SimInfo.Components.Sphere', num2str(targets(k)), '.componentType.charID;']);
                eval(['compNumber = m_MSMS.SimInfo.Components.Sphere', num2str(targets(k)), '.componentNumber.numID;']);
                eval(['nDOF       = m_MSMS.SimInfo.Components.Sphere', num2str(targets(k)), '.nDOF;']);
                
                % Builds the Sphere Feature
                % [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
                m_MSMS.udpPacket = [ m_MSMS.udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
                m_MSMS.udpPacket = [ m_MSMS.udpPacket; flipud( typecast( int16( compNumber ), 'int8' )' ) ]; % Joint NumID
                m_MSMS.udpPacket = [ m_MSMS.udpPacket; flipud( typecast( int16( nDOF ), 'int8' )' ) ]; % # Feature Number
                messLen = length( m_MSMS.udpPacket ) + 1;
                
                % Sets the [m or rad] chunk to 0.
                for n = 1:nDOF
                    m_MSMS.udpPacket = [ m_MSMS.udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero meters
                end % END FOR
                
                % Applies the position to the feature
                for m = 1:nDOF
                    m_MSMS.udpPacket( (messLen+4*(m-1)):(messLen+4*(m-1)+3) ) = flipud( typecast( single( m_MSMS.spherePos(targets(k),m) ), 'int8' )' );
                end % END FOR
                
            end % END FOR
            
            if singleton
                m_MSMS.SendUDP( );
            else
                % Debug
            end % END IF
            
        end % END FUNCTION
        
        % Colors the given Spheres
        function ColorSphere( m_MSMS, targets, intensity, singleton )
            
            % If no input for singleton is given, assume true
            if nargin < 4
                singleton = true;
            end % END IF
            
            % If function was called as a singleton, create and send the
            % UDP packet
            if singleton
                m_MSMS.CreateUDP( length(targets) * 1 );
                % DEBUG
            else
                % DEBUG out that func in not a singleton
            end % END IF
            
            % Parses Sphere info
            compType   = ' S';
            compNumber = [6, 10, 14, 18, 22, 26];
            nDOF       = 3;
            featureNum = 11;
            
            color = m_MSMS.colors .* intensity;
            
            for k = 1:length(targets)
                
                % Builds the Sphere Feature
                % [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
                m_MSMS.udpPacket = [ m_MSMS.udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
                m_MSMS.udpPacket = [ m_MSMS.udpPacket; flipud( typecast( int16( compNumber(targets(k)) ), 'int8' )' ) ]; % Joint NumID
                m_MSMS.udpPacket = [ m_MSMS.udpPacket; flipud( typecast( int16( featureNum ), 'int8' )' ) ]; % # Feature Number
                messLen = length( m_MSMS.udpPacket ) + 1;
                
                % Sets the color to 0|0|0.
                for n = 1:nDOF
                    m_MSMS.udpPacket = [ m_MSMS.udpPacket; int8(0) ]; % Zero meters
                end % END FOR
                
                % Applies the color to the feature
                for m = 1:nDOF
                    m_MSMS.udpPacket(messLen + m - 1) = flipud( typecast(  uint8(color(targets(k), m)) , 'uint8' )' );
                end % END FOR
                
            end % END FOR
            
            if singleton
                m_MSMS.SendUDP( );
            else
                % Debug
            end % END IF
            
        end % END FUNCTION
        
        % Adjusts the transparency of the given Spheres
        function HideSphere( m_MSMS, targets, transparency, singleton )
            
            % If no input for singleton is given, assume true
            if nargin < 4
                singleton = true;
            end % END IF
            
            % If function was called as a singleton, create and send the
            % UDP packet
            if singleton
                m_MSMS.CreateUDP( length(targets) * 1 );
                % DEBUG
            else
                % DEBUG out that func in not a singleton
            end % END IF
            
            
            compType   = ' S';
            compNumber = [6, 10, 14, 18, 22, 26];
            featureNum = 12;
            
            % Loops over the number of features
            for k = 1:length(targets)
                
                % Builds the packet.
                % [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
                m_MSMS.udpPacket = [ m_MSMS.udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
                m_MSMS.udpPacket = [ m_MSMS.udpPacket; flipud( typecast( int16( compNumber(targets(k)) ), 'int8' )' ) ]; % Joint NumID
                m_MSMS.udpPacket = [ m_MSMS.udpPacket; flipud( typecast( int16( featureNum ), 'int8' )' ) ]; % # Feature Number
                
                % Sets the trasnparency the given value.
                m_MSMS.udpPacket = [ m_MSMS.udpPacket; int8( transparency(k) ) ];
                
            end % END FOR
            
            if singleton
                m_MSMS.SendUDP( );
            else
                % Debug
            end % END IF
            
        end % END FUNCTION
        
        % Places the camera
        function PlaceCamera( m_MSMS, singleton )
            
            % If no input for singleton is given, assume true
            if nargin < 2
                singleton = true;
            end % END IF
            
            % If function was called as a singleton, create and send the
            % UDP packet
            if singleton
                m_MSMS.CreateUDP( 1 );
                % DEBUG
            else
                % DEBUG out that func in not a singleton
            end % END IF
            
            % Parses PlaceCamera info
            compType   = 'HT';
            compNumber = 0;
            nDOF       = 6;
            featureNum = 14;
            
            % Builds the WristRoll Feature
            % [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
            m_MSMS.udpPacket = [ m_MSMS.udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
            m_MSMS.udpPacket = [ m_MSMS.udpPacket; flipud( typecast( int16( compNumber ), 'int8' )' ) ]; % Joint NumID
            m_MSMS.udpPacket = [ m_MSMS.udpPacket; flipud( typecast( int16( featureNum ), 'int8' )' ) ]; % # Feature Number
            messLen = length( m_MSMS.udpPacket ) + 1;
            
            % Sets the [m or rad] chunk to 0.
            for n = 1:nDOF
                m_MSMS.udpPacket = [ m_MSMS.udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero all DOFs
            end % END FOR
            
            % XYZ-Pos: 0.5,0.09,0.45 [m]
            % XYZ-Rot: -30, 0.0, 0.0 [deg]
            % posRot = [ 0.50, (-30 / 360) * (2*pi) + (2*pi);...
            %            0.09, 0.0;...
            %            0.45, 0.0];
            
            % posRot = [ 0.71690, (-59.8471 / 360) * (2*pi) + (2*pi);...
            %            -0.0424, (34.722 / 360) * (2*pi);...
            %            0.02400, (-32.9075 / 360) * (2*pi) + (2*pi)];
            
            % posRot = [ 0.45000, (15 / 360) * (2*pi);...
            %            -0.2318, 0.0;...
            %            0.14250, 0.0];
            
            % posRot = [ 0.45000, (-75 / 360) * (2*pi) + (2*pi);...
            %            -0.0129, 0.0;...
            %            0.02100, 0.0];
            
            % posRot = [ 0.550, (-120 / 360) * (2*pi) + (2*pi);...
            %            0.091, 0.0;...
            %            -0.15, 0.0];
            
            % Applies the rotation to the feature
            for m = 1:nDOF
                m_MSMS.udpPacket( (messLen+4*(m-1)):(messLen+4*(m-1)+3) ) = flipud( typecast( single( m_MSMS.cameraPosRot(m) ), 'int8' )' );
            end % END FOR
            
            if singleton
                m_MSMS.SendUDP( );
            else
                % Debug
            end % END IF
            
        end % END FUNCTION
        
        % Move the given Joints 
        function MoveJoint( m_MSMS, joints, angle, targets, singleton )
            
            % If no input for singleton is given, assume true
            if nargin < 4
                singleton = true;
                targets = [1:6];
            elseif nargin < 5
                singleton = true;
            end% END IF
            
            % If function was called as a singleton, create and send the
            % UDP packet
            if singleton
                m_MSMS.CreateUDP( length(joints) * 1 );
                % DEBUG
            else
                % DEBUG out that func in not a singleton
            end % END IF
            
            for n = 1:length(targets)
            
                for k = 1:length(joints)
                    
                    % Parses ShoulderRoll info
                    compType   = eval(['m_MSMS.SimInfo.Components.', joints{k}, '.componentType.charID;']);
                    compNumber = eval(['m_MSMS.SimInfo.Components.', joints{k}, '.componentNumber.numID;']);
                    nDOF       = eval(['m_MSMS.SimInfo.Components.', joints{k}, '.nDOF;']);
                    
                    % Builds the ShoulderRoll Feature
                    % [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
                    m_MSMS.udpPacket = [ m_MSMS.udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
                    m_MSMS.udpPacket = [ m_MSMS.udpPacket; flipud( typecast( int16( compNumber ), 'int8' )' ) ]; % Joint NumID
                    m_MSMS.udpPacket = [ m_MSMS.udpPacket; flipud( typecast( int16( nDOF ), 'int8' )' ) ]; % # FeatureNum
                    messLen = length( m_MSMS.udpPacket ) + 1;
                    
                    % Sets the [m or rad] chunk to 0.
                    for n = 1:nDOF
                        m_MSMS.udpPacket = [ m_MSMS.udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % zero radians
                    end % END FOR
                    
                    % Applies the rotation to the feature
                    m_MSMS.udpPacket( (messLen):(messLen+3) ) = flipud( typecast( single( angle(n, k) ), 'int8' )' );
                    
                end % END FOR
            end % END FOR
            
            if singleton
                m_MSMS.SendUDP( );
            else
                % Debug
            end % END IF
            
        end % END FUNCTION
        
        % A proportional way to move fingers
        function MoveFingerProp( m_MSMS, inputStr )

            %%
            
            % Seperatees the input string into target/state pairs
            inputCell      = regexp( inputStr, m_MSMS.exprStr, 'Tokens' )';
            inputCell      = reshape([inputCell{:}],2,[])';
            inputCell(:,2) = cellfun(@(x) {str2double(x)},inputCell(:,2));
            
            sphereCell     = inputCell(~all(cellfun(@isempty, regexp(inputCell(:,1),'Sphere')), 2), :);
            targetCell     = inputCell( all(cellfun(@isempty, regexp(inputCell(:,1),'Sphere')), 2), :);
            
            m_MSMS.AngleFun( [targetCell{:,2} ])   
            
            %               DIP,PIP,Knuckle Joints      Hide and Color Spheres
            m_MSMS.CreateUDP( size(targetCell,1)*3 + size(sphereCell,1)*2 );
            
            
            if (m_MSMS.taskType ~= 1) && (m_MSMS.taskType ~= 2)
                m_MSMS.MoveJoint( {'ThumbDIP', 'ThumbKnuckle', 'ThumbSide'}, m_MSMS.jointAngle(1,1:3), [1], false );
            else
                m_MSMS.MoveJoint( {'ThumbDIP', 'ThumbKnuckle', 'ThumbMeta'}, m_MSMS.jointAngle(1,[1,2,4]), [1], false );
            end % END IF
            
            m_MSMS.MoveJoint( {'IndexDIP',  'IndexPIP',  'IndexKnuckle'},  m_MSMS.jointAngle(2, 1:3), [2], false );
            m_MSMS.MoveJoint( {'MiddleDIP', 'MiddlePIP', 'MiddleKnuckle'}, m_MSMS.jointAngle(3, 1:3), [3], false );
            m_MSMS.MoveJoint( {'RingDIP',   'RingPIP',   'RingKnuckle'},   m_MSMS.jointAngle(4, 1:3), [4], false );
            m_MSMS.MoveJoint( {'LittleDIP', 'LittlePIP', 'LittleKnuckle'}, m_MSMS.jointAngle(5, 1:3), [5], false );
            
            
            %%
            
            if m_MSMS.taskType == 3
                m_MSMS.HideSphere( [1,2,3,4,5], [sphereCell{:,2}].*100, false);
            end % END IF
            
            m_MSMS.SendUDP( );
                
            
        end % END FUNCTION
        
        % Calculates the angle for DIP/PIP/Knuckle joints from a given
        % proportion
        function AngleFun( m_MSMS, prop )
            
            % Flips the prop input to vertical
            if size(prop,1) == 1
                prop = prop(:);
            end
            
            % Repeats the matrix along the x-axis 4 times and multiplies it
            % by the differenece between angleMax and angleMin, then adds
            % angleMin on.
            m_MSMS.jointAngle = (repmat(prop, 1,4) .* (m_MSMS.jointAngleMax - m_MSMS.jointAngleMin)) + m_MSMS.jointAngleMin;
                
        end % END FUNCTION
        
        % Switches the taskType based on LabView input
        function TaskSwitch( m_MSMS, taskType )
            
            % Assigns the new taskType
            m_MSMS.taskType = taskType;
            
            % Sets up the scene with the chosen taskType
            m_MSMS.SetupScene( ) 
        end % END FUNCTION
          
        % END Working Bits
    end % END METHODS
    
end % END CLASS
