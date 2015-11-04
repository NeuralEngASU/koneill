%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:    RunDemo
% Desc:     Demos the capabilities of a MATLAB/MSMS pairing for running tasks
% Author:   Kevin O'Neill (Greger Lab)
% Date:     July 18, 2012
% 
% Use:      Open MSMS FingerTaskDemo
%           Setup and run animation
%           Type in "RunDemo" into MATLAB
%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RunDemo(  )

% Adds the correct folders to path
addpath('E:\MSMS\MSMS 2.1\bin\Matlab');
addpath('F:\CodeRepo\PatientCart\MSMS');
addpath('F:\CodeRepo\PatientCart\MSMS\Demos');
addpath('F:\CodeRepo\PatientCart\MSMS\PositionMat');
addpath('F:\CodeRepo\PatientCart\MSMS\Resources');

% Parses the XML file for the current demo
SimInfo = ParseXML('E:\MSMS\MSMS 2.1\Tutorials\FingerTaskDemo');

ParseRelPos( SimInfo );

SetupScene( SimInfo ) % Sets up the scene: Orientates camera, objects, arm, etc.
DemoTrialLoop( SimInfo ) % Loops through a demo of Finger Trials

end % END FUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - SetupScene - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetupScene( SimInfo )
% Sets up the scene: Places camera, objects, colors, etc.

% List of features to change during setup
featureList = {'ShoulderRoll'; 'ElbowPitch'; 'RotateWrist'; 'PlaceCamera';...
     'PlaceSphere1'; 'PlaceSphere2'; 'PlaceShpere3';...
     'ColorSphere1'; 'ColorSphere2'; 'ColorSphere3';...
     'HideSphere1';  'HideSphere2';  'HideSphere3'; ...
     'ZeroFingers'};

% Counts the number of features to change.
featureCount = length(featureList);

if strcmp(featureList{end},'ZeroFingers')
    featureCount = featureCount - 1;
end

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

% [~] = ModifyRelPos([1:6], 'ShoulderRoll', rotation);

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
% [~] = ModifyRelPos([1:6], 'ShoulderRoll', rotation);

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
% [~] = ModifyRelPos([1:6], 'ShoulderRoll', rotation);

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

%%%%% TMP %%%%%
% judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )

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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% - ZeroFingers - %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

pause(0.015)
MoveFingers( SimInfo, [1,2,3], 0);

end % END FUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% - DemoTrialLoop - %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DemoTrialLoop( SimInfo )

% Baseline data gathering. Patient should not move.
pause(5)

% loop for 150 trials
for k = 1:150
    
    % Generate random integers from 1 to 7. (Number of different trials)
    trialType = randint(1,1,[1,7]);
    
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
    end % END SWITCH
    
    % Makes the slected targets visable
    TargetVisable(targets,0);
    % Pause for 200 msec (~reaction time)
    pause(0.2)
    % Move fingers to targets
    MoveFingers(SimInfo,targets,1);
    % Hold position for 1 sec (Trial success determination
    pause(1)
    % Makes the targets invisable
    TargetVisable(targets,100);
    % pause for 200 msec (~reaction time)
    pause(0.2)
    % Move fingers to extended position
    MoveFingers(SimInfo,targets,0);
    % Pause for 1 sec for next trial
    pause(1)

end % END FOR

end % END FUNCTION


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ParseRelPos( SimInfo )


%%

pointStr = {'Thumb'; 'Index'; 'Middle'; 'Ring'; 'Little'; 'Palm'};

for k = 1:length(pointStr)-1
    eval(['relPos', pointStr{k}, ' = zeros(7, (4 + 1));']); % seven rows (xyz,theta,uvw) and columns for 3 finger joints and 1 POI.
end % END FOR

relPosPalm    = zeros(7, 1); % Only need offset from WristYaw
relPos        = zeros(7, (8 + 4)); % seven rows (xyz,theta,uvw) and columns for anchor point, 3 global rotations, and 8 joints.

relPosThumb(:,end)  = [ 0.0142; -0.0170; 0.0105; NaN; NaN; NaN; NaN];
relPosIndex(:,end)  = [-0.0095; -0.0140; 0.0052; NaN; NaN; NaN; NaN];
relPosMiddle(:,end) = [-0.0095; -0.0140; 0.0049; NaN; NaN; NaN; NaN];
relPosRing(:,end)   = [      0;       0;      0; NaN; NaN; NaN; NaN];
relPosLittle(:,end) = [      0;       0;      0; NaN; NaN; NaN; NaN];
relPosPalm(:,end)   = [      0;   -0.06;   0.04; NaN; NaN; NaN; NaN];

relPos(1:3,1) = [ 0.0900; 0.0501; 0.2599 ]; % Global Offset
relPos(4:7,1) = NaN;
relPos(5:7,2:4) = eye(3);

jointOrder = fileread('jointOrder.ini');

parseExpr{1} = '([A-Za-z0-9]+?)\s';
parseExpr{2} = '(Thumb[A-Za-z]+)';
parseExpr{3} = '(Index[A-Za-z]+)';
parseExpr{4} = '(Middle[A-Za-z]+)';
parseExpr{5} = '(Ring[A-Za-z]+)';
parseExpr{6} = '(Little[A-Za-z]+)';

parseArm    = regexp(jointOrder, parseExpr{1}, 'Tokens');
parseThumb  = regexp(jointOrder, parseExpr{2}, 'Tokens');
parseIndex  = regexp(jointOrder, parseExpr{3}, 'Tokens');
parseMiddle = regexp(jointOrder, parseExpr{4}, 'Tokens');
parseRing   = regexp(jointOrder, parseExpr{5}, 'Tokens');
parseLittle = regexp(jointOrder, parseExpr{6}, 'Tokens');

parseArm = parseArm(1:8);


for k = 1:8
    
    relPos(1:3, 4+k ) = SimInfo.Components.(char(parseArm{k})).abc';
    relPos(5:7, 4+k ) = SimInfo.Components.(char(parseArm{k})).uvw';
    
end % END FOR

relPos(4,5:10) = [0, (45/360)*(2*pi), 0, 0, (45/360)*(2*pi), (-30/360)*(2*pi)];

relPosMiddle(5:7,1) = [-1; 0; 0 ];

for k = 1:length(parseIndex)
    
    relPosThumb(1:3, k ) = SimInfo.Components.(char(parseThumb{k})).abc';
    relPosThumb(5:7, k ) = SimInfo.Components.(char(parseThumb{k})).uvw';
    
    relPosIndex(1:3, k ) = SimInfo.Components.(char(parseIndex{k})).abc';
    relPosIndex(5:7, k ) = SimInfo.Components.(char(parseIndex{k})).uvw';
    
    if k < length(parseIndex)
        relPosMiddle(1:3, 1+k ) = SimInfo.Components.(char(parseMiddle{k})).abc';
        relPosMiddle(5:7, 1+k ) = SimInfo.Components.(char(parseMiddle{k})).uvw';
    end 
    
    relPosRing(1:3, k ) = SimInfo.Components.(char(parseRing{k})).abc';
    relPosRing(5:7, k ) = SimInfo.Components.(char(parseRing{k})).uvw';

    relPosLittle(1:3, k ) = SimInfo.Components.(char(parseLittle{k})).abc';
    relPosLittle(5:7, k ) = SimInfo.Components.(char(parseLittle{k})).uvw';
    
end % END FOR


for k = 1:length(pointStr)
    eval(['relPos', pointStr{k}, ' = cat(2, relPos, relPos', pointStr{k},');']); % Concatenates relPos with each POI matrix.
end% END FOR

for k = 2:size(relPosIndex,2)
    
    if k < size(relPosPalm)
        relPosPalm(1:3,k) = relPosPalm(1:3,k-1) + relPosPalm(1:3,k);
    end
    
    relPosThumb(1:3,k)  = relPosThumb(1:3,k-1)  + relPosThumb(1:3,k);
    relPosIndex(1:3,k)  = relPosIndex(1:3,k-1)  + relPosIndex(1:3,k);
    relPosMiddle(1:3,k) = relPosMiddle(1:3,k-1) + relPosMiddle(1:3,k);
    relPosRing(1:3,k)   = relPosRing(1:3,k-1)   + relPosRing(1:3,k);
    relPosLittle(1:3,k) = relPosLittle(1:3,k-1) + relPosLittle(1:3,k);
    
end % END FOR

for k = 1:length(pointStr)
    eval(['relPos', pointStr{k}, ' = fliplr(relPos', pointStr{k},');']); % Flips the matricies from left to right.
end% END FOR

relPos = [];
relPos{1,1} = relPosThumb;
relPos{2,1} = relPosIndex;
relPos{3,1} = relPosMiddle;
relPos{4,1} = relPosRing;
relPos{5,1} = relPosLittle;
relPos{6,1} = relPosPalm;

% save('F:\CodeRepo\PatientCart\MSMS\PositionMat\relPosThumb.mat', 'relPosThumb');
% save('F:\CodeRepo\PatientCart\MSMS\PositionMat\relPosIndex.mat', 'relPosIndex');
% save('F:\CodeRepo\PatientCart\MSMS\PositionMat\relPosMiddle.mat','relPosMiddle');
% save('F:\CodeRepo\PatientCart\MSMS\PositionMat\relPosRing.mat',  'relPosRing');
% save('F:\CodeRepo\PatientCart\MSMS\PositionMat\relPosLittle.mat','relPosLittle');
% save('F:\CodeRepo\PatientCart\MSMS\PositionMat\relPosPalm.mat',  'relPosPalm');

save ('F:\CodeRepo\PatientCart\MSMS\PositionMat\relPos.mat', 'relPos');


end % END FUNCTION

function relPos = ModifyRelPos( target, jointList, angle )

    if ~iscell(jointList)
        jointList = cellstr(jointList);
    end

    load('F:\CodeRepo\PatientCart\MSMS\PositionMat\relPos.mat');
    
    % Joints considered for one target.
    jointLocalList  = {'POIOffset'; 'DIP'; 'PIP'; 'Knuckle'; 'Side'};
    % The thumb is really POIoffset, DIP, Knuckle, Yoke, Meta. But the
    % ordering is the same so it does not really matter.
    
    % Joints considered by all targets.
    jointSharedList = {'WristYaw';     'WristPitch';    'WristRoll';...
        'ElbowPitch';   'ShoulderRoll2'; 'ShoulderPitch';...
        'ShoulderRoll'; 'Shoulder1';     'GlobalXRot';...
        'GlobalYRot';   'GlobalZRot';    'GlobalOffset'};
    
    for j = 1:length(jointList)
        
        jointName = jointList{j};
        
        % TRUE if the joint in question is shared by all targets.
        sharedFlag = true;
        
        % Detects if the joint is local or global
        for k = 1:length(jointLocalList)
            if strcmp(jointName, char(jointLocalList{k}))
                sharedFlag = false;
            end % END IF
        end % END FOR
        
        % Based on the location of the joint, edit the necessary targets.
        if sharedFlag
            % Loops over jointSharedList
            for k = 1:length(jointSharedList)
                % Finds the joint specified
                if strcmp(jointName, char(jointSharedList{k}))
                    % Loops over all targets
                    for n = 1:size(relPos,1)
                        % Shared joints always start at 12 back from the end.
                        relPos{n,1}(4,end-(length(jointSharedList)-k)) = angle(j);
                    end % END FOR
                end % END IF
            end % END FOR
        else
            for k = 1:length(jointLocalList)
                if strcmp(jointName, char(jointLocalList{k}))
                    relPos{target,1}(4,k) = angle(j);
                end% END IF
            end % END FOR
        end% END IF
        
    end
    
    save('F:\CodeRepo\PatientCart\MSMS\PositionMat\relPos.mat', 'relPos');
    
    relPos = relPos{target,1};
    

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - MoveFingers - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MoveFingers( SimInfo, targets, state )

% Loads in the relative position matricies for each finger in order to
% calculate position.
load('F:\CodeRepo\PatientCart\MSMS\PositionMat\relPosIndex.mat');
load('F:\CodeRepo\PatientCart\MSMS\PositionMat\relPosMiddle.mat');
load('F:\CodeRepo\PatientCart\MSMS\PositionMat\relPosThumb.mat');

% compType and nDOF are the same for every joint
compTypeFinger = ' J';
nDOFFinger     = 1;

% Total number of loops
waitTime = 0.02;
maxLoop = 0.5/waitTime; %1/2 sec divided by smaller pauses

%%%%%%%%%%
% Matricies to contain the rotation steps
rot2(1,:) = linspace(0,(45/360)*(2*pi),maxLoop);
rot2(2,:) = linspace(0,(10/360)*(2*pi),maxLoop);

rot3(1,:) = linspace(0,(10/360)*(2*pi),maxLoop);
rot3(2,:) = linspace(0,(70/360)*(2*pi),maxLoop);
rot3(3,:) = linspace(0,(30/360)*(2*pi),maxLoop);

rotE(1,:) = linspace((-30/360)*(2*pi),(60/360)*(2*pi),maxLoop);
%%%%%%%%%%

% state = 1 means close fist
% state = 0 means open fist
% If fist will close, reverse the rotation-step-matricies
if ~state
    rot2 = fliplr(rot2);
    rot3 = fliplr(rot3);
    rotE = fliplr(rotE);
end

% Initializes featureCount
featureCount = 0;

% Counts the number of features to change. Thumb = 2, Index/Middle = 3;
for n = 1:length(targets)
    if targets(n) == 1
        featureCount = featureCount + 2 + 1; % 2 joints, 1 color
    else
        featureCount = featureCount + 3 + 1; % 3 joints, 1 color
    end
end

% Loops for length of time.
for k = 1:maxLoop
    
    tic1 = tic;
    
    % Initializes the udpPacket.
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
                    udpPacket = [ udpPacket; int8( compTypeFinger(1) ); int8( compTypeFinger(2) ) ]; % Joint CharID
                    udpPacket = [ udpPacket; flipud( typecast( int16( compNumber(n) ), 'int8' )' ) ]; % Joint NumID
                    udpPacket = [ udpPacket; flipud( typecast( int16( nDOFFinger ), 'int8' )' ) ]; % Feature Number
                    m(n) = length( udpPacket ) + 1;
                    
                    % Sets the [m or rad] PIP chunk to 0 radians or meters.
                    udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero radians
                    
                end
                
                % Inserts rotation values into the packet.
                udpPacket( (m(1)):(m(1)+3) ) = flipud( typecast( single( rot2(1,k) ), 'int8' )' );
                udpPacket( (m(2)):(m(2)+3) ) = flipud( typecast( single( rot2(2,k) ), 'int8' )' );
                
                %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% - Calculate Pos - %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  
                relPos = ModifyRelPos(finger, {'DIP'; 'PIP'; 'Knuckle'}, [ rot2(2,k), rot2(1,k), 0 ]);
                posPOI = CalculatePos(relPos);                 
%                 udpPacket = PlaceSphere(udpPacket, SimInfo, 1, posPOI);

                
                %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% - ColorSphere - %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Radius of Volume of Interest (in meters) ( 0.0075 m is
                % current radius of the sphere.
                radius = 1 * 0.0075;
                
                % Position of the sphere.
                posSphere = [ 0.5400, -0.1400, -0.0260 ];
                
                % Calculates the abs distance the point is from the center of the
                % sphere.
                difference = posSphere' - posPOI;
                distance = sqrt(sum(difference.^2));
                
                % Colors the sphere based on the distance and radius of VoI
                if distance > radius
                    color = 32;
                else
                    color = 127; %uint8((1-(distance/radius))*2);
                end
                
                udpPacket = ColorSphere( udpPacket, finger, color);
%                 
                
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
                    udpPacket = [ udpPacket; int8( compTypeFinger(1) ); int8( compTypeFinger(2) ) ]; % Joint CharID
                    udpPacket = [ udpPacket; flipud( typecast( int16( compNumber(n) ), 'int8' )' ) ]; % Joint NumID
                    udpPacket = [ udpPacket; flipud( typecast( int16( nDOFFinger ), 'int8' )' ) ]; % Feature Number
                    m(n) = length( udpPacket ) + 1;
                    
                    % Sets the [m or rad] PIP chunk to 0 radians or meters.
                    udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero radians
                    
                end
                
                % Inserts rotation values into the packet.
                udpPacket( (m(1)):(m(1)+3) ) = flipud( typecast( single( rot3(1,k) ), 'int8' )' );
                udpPacket( (m(2)):(m(2)+3) ) = flipud( typecast( single( rot3(2,k) ), 'int8' )' );
                udpPacket( (m(3)):(m(3)+3) ) = flipud( typecast( single( rot3(3,k) ), 'int8' )' );
                
                %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% - Calculate Pos - %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                relPos = ModifyRelPos(finger, {'DIP'; 'PIP'; 'Knuckle'}, [ rot3(3,k), rot3(2,k), rot3(1,k) ]);
                posPOI = CalculatePos(relPos);                 
                %%
                % Place Sphere
                
%                 udpPacket = PlaceSphere(udpPacket, SimInfo, 2, posPOI);
                
                %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% - ColorSphere - %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Radius of Volume of Interest (in meters) ( 0.0075 m is
                % current radius of the sphere.
                radius = 1 * 0.0075;
                
                % Position of the sphere.
                posSphere = [ 0.5750, -0.1435, -0.0110 ];
                
                % Calculates the abs distance the point is from the center of the
                % sphere.
                difference = posSphere' - posPOI;
                distance = sqrt(sum(difference.^2));
                
                % Colors the sphere based on the distance and radius of VoI
                if distance > radius
                    color = 32;
                else
                    color = 127; %uint8((1-(distance/radius))*2);
                end
                
                udpPacket = ColorSphere( udpPacket, finger, color);
                
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
                    udpPacket = [ udpPacket; int8( compTypeFinger(1) ); int8( compTypeFinger(2) ) ]; % Joint CharID
                    udpPacket = [ udpPacket; flipud( typecast( int16( compNumber(n) ), 'int8' )' ) ]; % Joint NumID
                    udpPacket = [ udpPacket; flipud( typecast( int16( nDOFFinger ), 'int8' )' ) ]; % Feature Number
                    m(n) = length( udpPacket ) + 1;
                    
                    % Sets the [m or rad] PIP chunk to 0 radians or meters.
                    udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero radians
                    
                end
                
                % Inserts rotation values into the packet.
                udpPacket( (m(1)):(m(1)+3) ) = flipud( typecast( single( rot3(1,k) ), 'int8' )' );
                udpPacket( (m(2)):(m(2)+3) ) = flipud( typecast( single( rot3(2,k) ), 'int8' )' );
                udpPacket( (m(3)):(m(3)+3) ) = flipud( typecast( single( rot3(3,k) ), 'int8' )' );
                
                %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% - Calculate Pos - %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                relPos = ModifyRelPos(finger, {'DIP'; 'PIP'; 'Knuckle'}, [ rot3(3,k), rot3(2,k), rot3(1,k) ]);
                posPOI = CalculatePos(relPos); 
                
%                 udpPacket = PlaceSphere(udpPacket, SimInfo, 3, posPOI);

                
                %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% - ColorSphere - %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Radius of Volume of Interest (in meters) ( 0.0075 m is
                % current radius of the sphere.
                radius = 1 * 0.0075;
                
                % Position of the sphere.
                posSphere = [ 0.5850, -0.1610, 0.0 ];
                
                % Calculates the abs distance the point is from the center of the
                % sphere.
                difference = posSphere' - posPOI;
                distance = sqrt(sum(difference.^2));
                
                % Colors the sphere based on the distance and radius of VoI
                if distance > radius
                    color = 32;
                else
                    color = 127; %uint8((1-(distance/radius))*2);
                end
                
                udpPacket = ColorSphere( udpPacket, finger, color);
                
            otherwise
        end % END SWITCH
    end % END FOR
    
    % Sends the packet
    judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )
    
    timeElapsed = toc(tic1);
    
    if timeElapsed < waitTime
        % Pauses for a small bit of time to allow the packet to register with
        % MSMS. MSMS reads packets at ~0.01 sec intervals
        pause(waitTime - timeElapsed)
    end
    
end % END FOR

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - CalculatePos - %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function posPOI = CalculatePos( relPos )

% Applies the new rotations to each joint in relPos
% relPos(4,2:4) = rot;

% The initial offset for the Point of Interest (POI) is the xyz-pos in
% the first column
offsetNew = relPos(1:3, 1);

% Iterate from the POI to the Anchor Point (AP). And rotate around
% the global axis rotation offset.
for n = 1: size(relPos,2) - 2
    
    % Pull out relevent info for current rotation
    theta = relPos(  4, n+1);
    abc   = relPos(1:3, n+1);
    xyz   = offsetNew;
    uvw   = relPos(5:7, n+1);
    
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
    offsetNew= xyzNew;
    
end % END FOR

% Add the POI position to the global offset
% posPOI = relPos(1:3,end) + offsetNew;
posPOI = offsetNew;


end % END FUNCTION


%%%%%%%
% Use for Debugging packets
%%%%%%%

% % % Initializes the udpPacket.
% udpPacket = []; % Empty packet
% udpPacket = [ udpPacket; int8(1) ];

%%%%% TMP %%%%%
% judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )