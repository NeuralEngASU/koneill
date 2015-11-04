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
    properties (SetAccess = private, GetAccess = private)
        % The 'm_MSMS.' before each variable is implicit.
        
        % Default params
        mode       = '0';
        taskType   = '0';
        errorLevel = '0';
        
        % Contains the structure for the parsed simulationSetup.xml and prostheses.xml
        simInfo   = [];
        
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
        
        % Rotations for each target
        rotMat = zeros(3,6,1);
        
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
        spherePos = zeros(6,3,4);
        
        % Position and rotation for the camera
        cameraPosRot = [ 0.50, (-30 / 360) * (2*pi) + (2*pi);...
                         0.09, 0.0;...
                         0.45, 0.0];
        %%%%%%%%%%% Keeping indent in check
        
        
        
    end % END PROPERTIES
    
    %% Enumeration of MSMS
    
    % Contains information for a few settings. Named variables makes it
    % easy for others to use this class.
    enumeration        
        % Male and Female models
        Male   (0)
        Female (1)
        
        % Left and Right arms
        Left   (0)
        Right  (1)
        
        % Quick and long calculation settings
        Light  (0)
        Full   (1)
        
        % Debug Log Levels
        Error   (0)
        Warning (1)
        Info    (2)
        
        % Points of Interest
        Thumb  (1)
        Index  (2)
        Middle (3)
        Ring   (4)
        Little (5)
        Palm   (6)
        
        % Spheres
        Sphere1 (1)
        Sphere2 (2)
        Sphere3 (3)
        Sphere4 (4)
        Sphere5 (5)
        Sphere6 (6)
        
        % Joint Information
        POIOffset     ( 1, 1)
        DIP           ( 2, 2)
        PIP           ( 3, 3)
        Knuckle       ( 4, 4)
        Side          ( 5, 5)
        WristYaw      ( 6, 0)
        WristPitch    ( 7, 0)
        WristRoll     ( 8, 0)
        ElbowPitch    ( 9, 0)
        ShoulderRoll2 (10, 0)
        ShoulderPitch (11, 0)
        ShoulderRoll  (12, 0)
        Shoulder1     (13, 0)
        GlobalXRot    (14, 6)
        GlobalYRot    (15, 7)
        GlobalZRot    (16, 8)
        GlobalOffset  (17, 9)  
    end % END ENUMERATION
    
    %% Methods of MSMS
    methods
        %% Object Creation/Initialization
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
        function init( m_MSMS, Params )
            
            m_MSMS.spherePos(:,:,1) = [ 0.5400, -0.1400, -0.0260;...
                                        0.5750, -0.1435, -0.0110;...
                                        0.5850, -0.1610,     0.0;...
                                           0.0,     0.0,     0.0;...
                                           0.0,     0.0,     0.0;...
                                           0.0,     0.0,     0.0];
            
            m_MSMS.spherePos(:,:,2) = [ 0.5400, -0.1400, -0.0260;...
                                        0.5750, -0.1435, -0.0110;...
                                        0.5850, -0.1610,     0.0;...
                                           0.0,     0.0,     0.0;...
                                           0.0,     0.0,     0.0;...
                                           0.0,     0.0,     0.0];
            
            
            % Creates and initializes the EventLog object.
            addpath('..\..\CoreFunctions\EventLog\');
            m_MSMS.m_eventLog = EventLog( );
            m_MSMS.m_eventLog.init( Params.Debug );
            
            % Parses Params
            m_MSMS.mode       = eval(['m_MSMS.', Params.mode, ';']);
            m_MSMS.taskType   = eval(['m_MSMS.', Params.taskType, ';']);
            m_MSMS.errorLevel = eval(['m_MSMS.', Params.errorLevel, ';']);
            m_MSMS.waitTime   = Params.waitTime;
            m_MSMS.loopTime   = Params.loopTime;
            m_MSMS.stateSteps = m_MSMS.loopTime/m_MSMS.waitTime;
            
            % Calculates the rotation value matrix
            m_MSMS.CalcRotMat();
            
            % Parses the simulationSetup.xml and prostheses.xml files from
            % the MSMS Scenario
            m_MSMS.SimInfo = ParseXML( Params.scenarioPath );
            
            % Parses and organizes data to form the relPos matrix.
            m_MSMS.ParseRelPos( );
            
            % Sets up the scene with the specified taskType
            m_MSMS.SetupScene( );
            
        end % END INIT
        
        % END Object Creation/Initialization
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
        
        % Gets the current 'relative Position' cell (relPos)
        function relPos = getRelPos( m_MSMS )
            relPos = m_MSMS.relPos;
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
        function SetupScene( m_MSMS )
            
            switch m_MSMS.taskType
                case 'Finger'
                    % Perhaps a 2 column cell, col1 for feature name, col2
                    % for feature count.
                    m_MSMS.featureList = {'ShoulderRoll', 1;   
                                            'ElbowPitch', 1;
                                           'RotateWrist', 1;
                                           'PlaceCamera', 1;
                                           'PlaceSphere', 3; %6;
                                           'ColorSphere', 3; %6;
                                            'HideSphere', 3; %6; 
                                           'ZeroFingers', 8}; %14};
                otherwise                    
            end % END SWITCH
            
            % Creates a UDP packet with the selected featureList
            m_MSMS.CreateUDP( sum(m_MSMS.featureList{:,2}) );
            
%%%%%%%%%%%%% Will need to think of a way to parse the featureList

            m_MSMS.MoveJoint( {'ShoulderRoll', 'ElbowPitch', 'WristRoll' },...
                [(45/360)*(2*pi), (45/360)*(2*pi), (-30/360)*(2*pi) + 2*pi],...
                [1:6], false) 
            
            m_MSMS.PlaceCamera( );
            m_MSMS.PlaceSphere( [1,2,3], false );
            m_MSMS.ColorSphere( [1,2,3], 127, false );
            m_MSMS.HideSphere(  [1,2,3], 100, false );
            m_MSMS.MoveFinger( );
            
%%%%%%%%%%%%%

            m_MSMS.SendUDP;
            
        end % END FUNCTION
        
        % Modifys selected joints from the Relative Position matrix by the
        % given angle.
        function ModifyRelPos( m_MSMS, POI, jointList, angleList )
            % Iterates over the joints in the jointList
            for k = 1:length(jointList)
                % Depending on the calculation mode, calculate for the full
                % arm or from the wrist to the POIs
                switch m_MSMS.mode
                    case 0
                        m_MSMS.relPos{POI}(4,eval(['m_MSMS.', jointList{n}, '(1,2)'])) = angleList(n);
                    case 1
                        m_MSMS.relPos{POI}(4,eval(['m_MSMS.', jointList{n}, '(1,1)'])) = angleList(n);
                    otherwise
                end % END SWITCH
            end % END FOR
        end % END FUNCTION  
        
        % Calculates the position of the POIs in three-space
        function CalculatePos( m_MSMS, POI)
                   
            % Iterates over the number of Points of Interest
            % Look into making this use better matricx multiplication
            for k = 1:length(POI)
                
                relPosPOI = m_MSMS.relPos{k};
                
                % The initial offset for the Point of Interest (POI) is the xyz-pos in
                % the first column
                offsetNew = relPosPOI(1:3, 1);
                
                % Iterate from the POI to the Anchor Point (AP). And rotate around
                % the global axis rotation offset.
                for n = 1: size(relPosPOI,2) - 2
                    
                    % Pull out relevent info for current rotation
                    theta = relPosPOI(  4, n+1);
                    abc   = relPosPOI(1:3, n+1);
                    xyz   = offsetNew;
                    uvw   = relPosPOI(5:7, n+1);
                    
                    % Parse info into seperate variables (It makes it easier to
                    %                                     look at)
                    a = abc(1);
                    b = abc(2);
                    c = abc(3);
                    
                    x = xyz(1);
                    y = xyz(2);
                    z = xyz(3);
                    
                    u = uvw(1);
                    v = uvw(2);
                    w = uvw(3);
                    
                    % Rotates the point <x,y,z> around the vector <u,v,w>, which passes
                    % through the point <a,b,c>, by angle theta in radians.
                    xyzNew = [(a*(v^2 + w^2) - u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(theta)) + x*cos(theta) + (-c*v+b*w-w*y+v*z)*sin(theta);...
                              (b*(u^2 + w^2) - v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(theta)) + y*cos(theta) + ( c*u-a*w+w*x-u*z)*sin(theta);...
                              (c*(u^2 + v^2) - w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(theta)) + z*cos(theta) + (-b*u+a*v-v*x+u*y)*sin(theta)];
                    
                    % New offset of POI is now xyzNew
                    offsetNew = xyzNew;
                    
                end % END FOR
                
                % Add the POI position to the global offset
                % posPOI = relPos(1:3,end) + offsetNew;
                m_MSMS.posPOI(k) = offsetNew;
            end
            
            
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
            if nargin == 2
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
                    m_MSMS.udpPacket( (messLen+4*(m-1)):(messLen+4*(m-1)+3) ) = flipud( typecast( single( m_MSMS.sphereLoc(targets(k),m) ), 'int8' )' );
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
            if nargin == 2
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
            compNumber = [98, 2, 106];
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
            if nargin == 2
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
            compNumber = [98,2,106];
            featureNum = 12;
            
            % Loops over the number of features
            for k = 1:length(targets)
                
                % Builds the packet.
                % [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
                m_MSMS.udpPacket = [ m_MSMS.udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
                m_MSMS.udpPacket = [ m_MSMS.udpPacket; flipud( typecast( int16( compNumber(targets(k)) ), 'int8' )' ) ]; % Joint NumID
                m_MSMS.udpPacket = [ m_MSMS.udpPacket; flipud( typecast( int16( featureNum ), 'int8' )' ) ]; % # Feature Number
                
                % Sets the trasnparency the given value.
                m_MSMS.udpPacket = [ m_MSMS.udpPacket; int8( transparency ) ];
                
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
            if nargin == 2
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
                m_MSMS.udpPacket = [ m_MSMS.udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; %#ok<AGROW> % zero radians
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
                m_MSMS.udpPacket( (messLen+4*(m-1)):(messLen+4*(m-1)+3) ) = flipud( typecast( single( posRot(m) ), 'int8' )' );
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
            if nargin <= 3
                singleton = true;
            elseif nargin <= 2
                singleton = true;
                targets = [1:6];
            end% END IF
            
            % If function was called as a singleton, create and send the
            % UDP packet
            if singleton
                m_MSMS.CreateUDP( length(joints) * 1 );
                % DEBUG
            else
                % DEBUG out that func in not a singleton
            end % END IF
            
            for k = 1:length(joints)
                
                % Parses ShoulderRoll info
                compType   = eval(['SimInfo.Components.', joints{k}, '.componentType.charID;']);
                compNumber = eval(['SimInfo.Components.', joints{k}, '.componentNumber.numID;']);
                nDOF       = eval(['SimInfo.Components.', joints{k}, '.nDOF;']);
                
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
                
                % Modifies the relPos matrix
                m_MSMS.ModifyRelPos( targets, joints{k}, angle(k) );
                
                % Applies the rotation to the feature
                m_MSMS.udpPacket( (messLen):(messLen+3) ) = flipud( typecast( single( angle(k) ), 'int8' )' );
                
            end % END FOR
            
            if singleton
                m_MSMS.SendUDP( );
            else
                % Debug
            end % END IF
            
        end % END FUNCTION
        
        % A state based way to move fingers
        function MoveFingerState( m_MSMS, inputStr )
            

        end % END FUNCTION
        
        % A proportional way to move fingers
        function MoveFingerProp( m_MSMS, targets, prop )

            % Seperatees the input string into target/state pairs
            inputCell      = regexp( inputStr, m_MSMS.exprStr, 'Tokens' )';
            inputCell      = reshape([inputCell{:}],2,[])';
            inputCell(:,2) = cellfun(@(x) {str2double(x)},inputCell(:,2));
            
            sphereCell     = inputCell(~all(cellfun(@isempty, regexp(inputCell(:,1),'Sphere')), 2), :);
            targetCell     = inputCell( all(cellfun(@isempty, regexp(inputCell(:,1),'Sphere')), 2), :);
            
            for k = 1:m_MSMS.stateSteps
                
                tic(clock1)
                
                %               DIP,PIP,Knuckle Joints      Hide and Color Spheres
                m_MSMS.CreateUDP( size(targetCell,1)*3 + size(sphereCell,1)*2 );
                
                m_MSMS.MoveJoint( {'DIP', 'PIP', 'Knuckle'}, [cellfun(@(x) eval(['m_MSMS.', x]))], [targetCell{:,2}], false );
                m_MSMS.HideSphere( [cellfun(@(x) eval(['m_MSMS.', x]))], [sphereCell{:,2}*100], false);
                
                % Determine success for all TRUE spheres
                % If success add intensity to sphere 
                % If failure subtract intensity from sphere
                
                m_MSMS.SendUDP( );
                
                packetTime = toc(clock1);
                
                if packetTime < m_MSMS.waitTime
                    pause( m_MSMS.waitTime - packetTime )
                end % END IF
                
            end % END FOR
            
            
        end % END FUNCTION
        
        % Detects if the target was sucessfully hit
        function DetectSucess( m_MSMS, target )
            
        end % END FUNCTION
          
        % Calculates the rotational value matrix
        function CalcRotMat( m_MSMS )
            
            % Pre allocate matrix
            m_MSMS.rotMat = zeros(3,6,m_MSMS.stateSteps);
            
            % Thumb
            m_MSMS.rotMat(1,1,:) = linspace(0,(45/360)*(2*pi),m_MSMS.stateSteps);
            m_MSMS.rotMat(2,1,:) = linspace(0,(10/360)*(2*pi),m_MSMS.stateSteps);
   
            for k = 2:4
                m_MSMS.rotMat(1,k,:) = linspace(0,(10/360)*(2*pi),m_MSMS.stateSteps);
                m_MSMS.rotMat(2,k,:) = linspace(0,(70/360)*(2*pi),m_MSMS.stateSteps);
                m_MSMS.rotMat(3,k,:) = linspace(0,(30/360)*(2*pi),m_MSMS.stateSteps);
            end % END FOR
            
        end % END FUNCTION
        % END Working Bits
    end % END METHODS
    
end % END CLASS
