%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:    RunDemo
% Desc:     Demos the capabilities of a MATLAB/MSMS pairing for running tasks
% Author:   Kevin O'Neill (Greger Lab)
% Date:     July 24, 2012
% 
% Use:      Open MSMS FingerTaskDemo
%           Setup and run animation
%           Add the D:\CodeRepo\PatientCart\MATLAB\MSMS and E:\MSMS folders to MATLAB's path
%           Type in "RunDemo" into MATLAB
%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% - DemoTrialLoop - %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DemoRTFRD( )

addpath('E:\MSMS\MSMS 2.1\bin\Matlab');
addpath('F:\CodeRepo\PatientCart\MSMS');
addpath('F:\CodeRepo\PatientCart\MSMS\Demos');
addpath('F:\CodeRepo\PatientCart\MSMS\PositionMat');
addpath('F:\CodeRepo\PatientCart\MSMS\Resources');

% Parses the XML file for the current demo
SimInfo = ParseXML('E:\MSMS\MSMS 2.1\Tutorials\FingerTaskDemo');

if ~isempty(instrfind('Type','udp')); fclose(instrfind('Type','udp')); end

u = udp('127.0.0.1',4201,'localhost','127.0.0.1','localport',4200);
fopen(u);
if u.BytesAvailable; flushinput(u); flushouput(u); end %clear buffers

%[1-7];0/1 | [1-7]; 2/3 (sphere off/sphere on)


while true
    
    paramStr = fscanf(u);
    
    if ~isempty(paramStr)
        paramCell = regexp(paramStr, ';', 'Split');
        
        trialType = round(str2double(paramCell{1,1}));
        state = round(str2double(paramCell{1,2}));
        
        switch trialType
            case 1
                targets = [1];
            case 2
                targets = [2];
            case 3
                targets = [3];
            case 4
                targets = [1,2];
            case 5
                targets = [1,3];
            case 6
                targets = [2,3];
            case 7
                targets = [1,2,3];
            otherwise
                SetupScene(SimInfo);
                targets = [];
        end % END SWITCH
        
        if (state == 0) || (state == 1)
            MoveFingers(SimInfo,targets,state);
        elseif (state == 2) || (state == 3)
            TargetVisable(targets, abs(state-3)*100);
        else
        end
        
    end % END IF
    
end % END WHILE

end % END FUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - SetupScene - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetupScene( SimInfo )
% Sets up the scene: Places camera, objects, colors, etc.

% List of features to change during setup
featureList = {'ShoulderRoll'; 'ElbowPitch'; 'RotateWrist'; 'PlaceCamera';...
     'PlaceSphere1'; 'PlaceSphere2'; 'PlaceShpere3';...
     'ColorSphere1'; 'ColorSphere2'; 'ColorSphere3';...
     'HideSphere1';  'HideSphere2';  'HideSphere3'};

% Counts the number of features to change.
featureCount = length(featureList);

% Palm pos = 0.4503, -0.1863, -0.0258 m
% Palm vector = 0,0,-1
% Palm rot = 0, 330, 90 deg


% Initializes the udpPacket.
udpPacket = []; % Empty packet
udpPacket = [ udpPacket; int8(featureCount) ];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% - ShoulderRoll - %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parses ShoulderRoll info
compType   = SimInfo.Components.ShoulderRoll.componentType.charID;
compNumber = SimInfo.Components.ShoulderRoll.componentNumber.numID;
nDOF       = SimInfo.Components.ShoulderRoll.nDOF;

% Builds the ShoulderRoll Feature
% [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
udpPacket = [ udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
udpPacket = [ udpPacket; flipud( typecast( int16( compNumber ), 'int8' )' ) ]; % Joint NumID
udpPacket = [ udpPacket; flipud( typecast( int16( nDOF ), 'int8' )' ) ]; % # FeatureNum
messLen = length( udpPacket ) + 1;

% Sets the [m or rad] chunk to 0.
for n = 1:nDOF
    udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; %#ok<AGROW> % zero radians
end % END FOR

% Rotation: 45 deg
rotation = 0.125*(2*pi);

% Applies the rotation to the feature
udpPacket( (messLen):(messLen+3) ) = flipud( typecast( single( rotation ), 'int8' )' );

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% - ElbowPitch - %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parses ElbowPitch info
compType   = SimInfo.Components.ElbowPitch.componentType.charID;
compNumber = SimInfo.Components.ElbowPitch.componentNumber.numID;
nDOF       = SimInfo.Components.ElbowPitch.nDOF;

% Builds the ElbowPitch Feature
% [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
udpPacket = [ udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
udpPacket = [ udpPacket; flipud( typecast( int16( compNumber ), 'int8' )' ) ]; % Joint NumID
udpPacket = [ udpPacket; flipud( typecast( int16( nDOF ), 'int8' )' ) ]; % # FeatureNum
messLen = length( udpPacket ) + 1;

% Sets the [m or rad] chunk to 0.
for n = 1:nDOF
    udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; %#ok<AGROW> % zero radians
end % END FOR

% Rotation: 45 deg
rotation = 0.125*(2*pi);

% Applies the rotation to the feature
udpPacket( (messLen):(messLen+3) ) = flipud( typecast( single( rotation ), 'int8' )' );

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% - WristRoll - %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parses WristRoll info
compType   = SimInfo.Components.WristRoll.componentType.charID;
compNumber = SimInfo.Components.WristRoll.componentNumber.numID;
nDOF       = SimInfo.Components.WristRoll.nDOF;

% Builds the WristRoll Feature
% [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
udpPacket = [ udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
udpPacket = [ udpPacket; flipud( typecast( int16( compNumber ), 'int8' )' ) ]; % Joint NumID
udpPacket = [ udpPacket; flipud( typecast( int16( nDOF ), 'int8' )' ) ]; % # FeatureNum
messLen = length( udpPacket ) + 1;

% Sets the [m or rad] chunk to 0.
for n = 1:nDOF
    udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; %#ok<AGROW> % zero radians
end % END FOR

% Rotation: -30 deg
rotation = (-30 / 360) * (2*pi) + (2*pi);

% Applies the rotation to the feature
udpPacket( (messLen):(messLen+3) ) = flipud( typecast( single( rotation ), 'int8' )' );

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% - PlaceCamera - %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % Initializes the udpPacket.
% udpPacket = []; % Empty packet
% udpPacket = [ udpPacket; int8(1) ];

% Parses PlaceCamera info
compType   = 'HT';
compNumber = 0;
nDOF       = 6;
featureNum = 14;

% Builds the WristRoll Feature
% [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
udpPacket = [ udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
udpPacket = [ udpPacket; flipud( typecast( int16( compNumber ), 'int8' )' ) ]; % Joint NumID
udpPacket = [ udpPacket; flipud( typecast( int16( featureNum ), 'int8' )' ) ]; % # Feature Number
messLen = length( udpPacket ) + 1;

% Sets the [m or rad] chunk to 0.
for n = 1:nDOF
    udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; %#ok<AGROW> % zero radians
end % END FOR

% XYZ-Pos: 0.5,0.09,0.45 [m]
% XYZ-Rot: -30, 0.0, 0.0 [deg]
posRot = [ 0.50, (-30 / 360) * (2*pi) + (2*pi);...
           0.09, 0.0;...
           0.45, 0.0];
       
% posRot = [ 0.7169, (-59.8471 / 360) * (2*pi) + (2*pi);...
%            -0.0424, (34.722 / 360) * (2*pi);...
%            0.024, (-32.9075 / 360) * (2*pi) + (2*pi)];
       
% posRot = [ 0.45, (15 / 360) * (2*pi);...
%            -0.2318, 0.0;...
%            0.1425, 0.0];
       
% posRot = [ 0.45, (-75 / 360) * (2*pi) + (2*pi);...
%           -0.0129, 0.0;...
%            0.021, 0.0];
       
% posRot = [ 0.55, (-120 / 360) * (2*pi) + (2*pi);...
%            0.091, 0.0;...
%            -0.15, 0.0];

% Applies the rotation to the feature
for m = 1:nDOF
    udpPacket( (messLen+4*(m-1)):(messLen+4*(m-1)+3) ) = flipud( typecast( single( posRot(m) ), 'int8' )' );
end % END FOR


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% - PlaceSphere - %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
udpPacket = PlaceSphere(udpPacket, SimInfo, [1,2,3]);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% - ColorSphere - %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
udpPacket = ColorSphere(udpPacket, [1,2,3], 127);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% - HideSphere - %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
udpPacket = HideSphere(udpPacket, [1,2,3], 100);

%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%% - SendUDP - %%%%%
%%%%%%%%%%%%%%%%%%%%%%%
judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )


end % END FUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - MoveFingers - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MoveFingers( SimInfo, targets, state )

if isempty(targets)
    return
end

% compType and nDOF are the same for every joint
compType = ' J';
nDOF     = 1;

% Total number of loops
maxLoop = 0.5/0.015; %1/2 sec divided by smaller pauses

%%%%%%%%%%
% Matricies to contain the rotation steps
rot2(1,:) = linspace(0,(45/360)*(2*pi),maxLoop);
rot2(2,:) = linspace(0,(10/360)*(2*pi),maxLoop);

rot3(1,:) = linspace(0,(10/360)*(2*pi),maxLoop);
rot3(2,:) = linspace(0,(70/360)*(2*pi),maxLoop);
rot3(3,:) = linspace(0,(30/360)*(2*pi),maxLoop);
%%%%%%%%%%

% state = 1 means close fist
% state = 0 means open fist
% If fist will close, reverse the rotation-step-matricies
if ~state
    rot2 = fliplr(rot2);
    rot3 = fliplr(rot3);
end

% Initializes featureCount
featureCount = 0;

% Counts the number of features to change. Thumb = 2, Index/Middle = 3;
for n = 1:length(targets)
    if targets(n) == 1
        featureCount = featureCount+2;
    else
        featureCount = featureCount +3;
    end
end

% Loops for length of time.
for k = 1:maxLoop
    
% % Initializes the udpPacket.
udpPacket = []; % Empty packet
udpPacket = [ udpPacket; int8(featureCount) ];

    % Loops over the targets/fingers
    for j = 1:length(targets)
        
        % Builds the packet with by successivly adding fingers.
        finger = targets(j);
        switch finger
            case 1
                % Parses the ThumbKnuckle information
                compNumber(1) = SimInfo.Components.ThumbKnuckle.componentNumber.numID;
                % Parses the ThumbDIP information
                compNumber(2) = SimInfo.Components.ThumbDIP.componentNumber.numID;
                
                % Loops over each joint
                for n = 1:2
                    
                    % Builds the packet
                    % [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
                    udpPacket = [ udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
                    udpPacket = [ udpPacket; flipud( typecast( int16( compNumber(n) ), 'int8' )' ) ]; % Joint NumID
                    udpPacket = [ udpPacket; flipud( typecast( int16( nDOF ), 'int8' )' ) ]; % Feature Number
                    m(n) = length( udpPacket ) + 1;
                    
                    % Sets the [m or rad] PIP chunk to 0 radians or meters.
                    udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero radians
                    
                end
                
                % Inserts rotation values into the packet.
                udpPacket( (m(1)):(m(1)+3) ) = flipud( typecast( single( rot2(1,k) ), 'int8' )' );
                udpPacket( (m(2)):(m(2)+3) ) = flipud( typecast( single( rot2(2,k) ), 'int8' )' );
            case 2            
                % Parses the IndexKnuckle information
                compNumber(1) = SimInfo.Components.IndexKnuckle.componentNumber.numID;
                % Parses the IndexPIP information
                compNumber(2) = SimInfo.Components.IndexPIP.componentNumber.numID;
                % Parses the IndexDIP information
                compNumber(3) = SimInfo.Components.IndexDIP.componentNumber.numID;
                
                % Loops over each joint
                for n = 1:3
                    
                    % Builds the packet
                    % [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
                    udpPacket = [ udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
                    udpPacket = [ udpPacket; flipud( typecast( int16( compNumber(n) ), 'int8' )' ) ]; % Joint NumID
                    udpPacket = [ udpPacket; flipud( typecast( int16( nDOF ), 'int8' )' ) ]; % Feature Number
                    m(n) = length( udpPacket ) + 1;
                    
                    % Sets the [m or rad] PIP chunk to 0 radians or meters.
                    udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero radians
                    
                end
                
                % Inserts rotation values into the packet.
                udpPacket( (m(1)):(m(1)+3) ) = flipud( typecast( single( rot3(1,k) ), 'int8' )' );
                udpPacket( (m(2)):(m(2)+3) ) = flipud( typecast( single( rot3(2,k) ), 'int8' )' );
                udpPacket( (m(3)):(m(3)+3) ) = flipud( typecast( single( rot3(3,k) ), 'int8' )' );
            case 3
                % Parses the MiddleKnuckle information
                compNumber(1) = SimInfo.Components.MiddleKnuckle.componentNumber.numID;
                % Parses the MiddlePIP information
                compNumber(2) = SimInfo.Components.MiddlePIP.componentNumber.numID;
                % Parses the MiddleDIP information
                compNumber(3) = SimInfo.Components.MiddleDIP.componentNumber.numID;
                
                % Loops over each joint
                for n = 1:3
                    
                    % Builds the packet
                    % [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
                    udpPacket = [ udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
                    udpPacket = [ udpPacket; flipud( typecast( int16( compNumber(n) ), 'int8' )' ) ]; % Joint NumID
                    udpPacket = [ udpPacket; flipud( typecast( int16( nDOF ), 'int8' )' ) ]; % Feature Number
                    m(n) = length( udpPacket ) + 1;
                    
                    % Sets the [m or rad] PIP chunk to 0 radians or meters.
                    udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero radians
                    
                end
                
                % Inserts rotation values into the packet.
                udpPacket( (m(1)):(m(1)+3) ) = flipud( typecast( single( rot3(1,k) ), 'int8' )' );
                udpPacket( (m(2)):(m(2)+3) ) = flipud( typecast( single( rot3(2,k) ), 'int8' )' );
                udpPacket( (m(3)):(m(3)+3) ) = flipud( typecast( single( rot3(3,k) ), 'int8' )' );
            otherwise
        end % END SWITCH
    end % END FOR
    
    % Sends the packet
    judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )
    
    % Pauses for a small bit of time to allow the packet to register with
    % MSMS. MSMS reads packets at ~0.01 sec intervals
    pause(0.015)
    
end % END FOR

end % END FUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% - TargetVisable - %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TargetVisable( targets, transparency )

% targets = [1,2,3];
% targets = [2,3];
% targets = [2];
% transparency = 0;
% transparency = 50;
% transparency = 100;

% Number of features to change.
featureCount = length(targets);

% compType and featureNum are the same for eavh feature (in this case).
% Sphere contains the compNumber for each sphere (Sphere1,2,3...)
compType   = ' S';
sphere = [98,2,106];
featureNum = 12;

% % Initializes the udpPacket.
udpPacket = []; % Empty packet
udpPacket = [ udpPacket; int8(featureCount) ];

% Loops over the number of features
for k = 1:featureCount
    
    % Extracts the compNumber from sphere
    compNumber = sphere(targets(k));

    % Builds the packet.
    % [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
    udpPacket = [ udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
    udpPacket = [ udpPacket; flipud( typecast( int16( compNumber ), 'int8' )' ) ]; % Joint NumID
    udpPacket = [ udpPacket; flipud( typecast( int16( featureNum ), 'int8' )' ) ]; % # Feature Number
    
    % Sets the trasnparency the given value.
    udpPacket = [ udpPacket; int8( transparency ) ];
    
end % END FOR

% Sends the packet.
judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )

end % END FUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - PlaceSphere - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function udpPacket = PlaceSphere( udpPacket, SimInfo, targets, location ) 
% Places the spheres given.

if nargin < 4
    location = [ 0.5400, -0.1400, -0.0260;...
                 0.5750, -0.1435, -0.0110;...
                 0.5850, -0.1610, 0.0 ];
else
    locationTmp = zeros(3);
    locationTmp(targets,:) = location(:)';
    location = locationTmp;
end
    

for k=1:length(targets)  
    
    
    % Parses Sphere info
    eval(['compType   = SimInfo.Components.Sphere', num2str(targets(k)), '.componentType.charID;']);
    eval(['compNumber = SimInfo.Components.Sphere', num2str(targets(k)), '.componentNumber.numID;']);
    eval(['nDOF       = SimInfo.Components.Sphere', num2str(targets(k)), '.nDOF;']);
    
    % Builds the Sphere Feature
    % [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
    udpPacket = [ udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
    udpPacket = [ udpPacket; flipud( typecast( int16( compNumber ), 'int8' )' ) ]; % Joint NumID
    udpPacket = [ udpPacket; flipud( typecast( int16( nDOF ), 'int8' )' ) ]; % # Feature Number
    messLen = length( udpPacket ) + 1;
    
    % Sets the [m or rad] chunk to 0.
    for n = 1:nDOF
        udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero meters
    end % END FOR
    
    % Applies the position to the feature
    for m = 1:nDOF
        udpPacket( (messLen+4*(m-1)):(messLen+4*(m-1)+3) ) = flipud( typecast( single( location(targets(k),m) ), 'int8' )' );
    end % END FOR
    
end % END FOR

end % END FUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - ColorSphere - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function udpPacket = ColorSphere( udpPacket, targets, color ) 
% Colors the spheres given.

% Parses Sphere info
compType   = ' S';
compNumber = [98, 2, 106];
nDOF       = 3;
featureNum = 11;

rgb = eye(3) .* color;

for k = 1:length(targets)
    
    % Builds the Sphere2 Feature
    % [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
    udpPacket = [ udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
    udpPacket = [ udpPacket; flipud( typecast( int16( compNumber(targets(k)) ), 'int8' )' ) ]; % Joint NumID
    udpPacket = [ udpPacket; flipud( typecast( int16( featureNum ), 'int8' )' ) ]; % # Feature Number
    messLen = length( udpPacket ) + 1;
    
    % Sets the color to 0|0|0.
    for n = 1:nDOF
        udpPacket = [ udpPacket; int8(0) ]; % Zero meters
    end % END FOR
    
    % Applies the color to the feature
    for m = 1:nDOF
        udpPacket(messLen + m - 1) = flipud( typecast(  uint8(rgb(targets(k), m)) , 'uint8' )' );
    end % END FOR

end % END FOR

end % END FUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - HideSphere - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function udpPacket = HideSphere( udpPacket, targets, transparency)

compType   = ' S';
compNumber = [98,2,106];
featureNum = 12;

% Loops over the number of features
for k = 1:length(targets)

    % Builds the packet.
    % [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
    udpPacket = [ udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
    udpPacket = [ udpPacket; flipud( typecast( int16( compNumber(targets(k)) ), 'int8' )' ) ]; % Joint NumID
    udpPacket = [ udpPacket; flipud( typecast( int16( featureNum ), 'int8' )' ) ]; % # Feature Number
    
    % Sets the trasnparency the given value.
    udpPacket = [ udpPacket; int8( transparency ) ];
    
end % END FOR

end % END FUNCTION


%%%%%%%
% Use for Debugging packets
%%%%%%%

% % % Initializes the udpPacket.
% udpPacket = []; % Empty packet
% udpPacket = [ udpPacket; int8(1) ];

%%%%% TMP %%%%%
% judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )