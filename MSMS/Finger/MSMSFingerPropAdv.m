%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:    MSMSFingerPropAdv
% Desc:     Sends the signals to MSMS to proportionally control finger joints
% Author:   Kevin O'Neill (Greger Lab)
% Date:     July 10, 2012
% 
% Use:      MSMSFingerProp( SimulationInfo, prop ) - Will move fingers proportionally from an open to closed state
%           MSMSFingerProp( SimulationInfo, [5,5,NaN,NaN,5] ) - Will fully
%           close the Index, Middle, and Thumb.
%           MSMSFingerProp( SimulationInfo, 1 )   - Bends the ElbowPitch joint by 90 deg
%
%           "prop" must be in order of [Index, Middle, Ring, Little, Thumb]
%           (due to seqNum) and must be either 5 or 1 in length.
%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function MSMSFingerPropAdv( SimulationInfo, prop )

% If prop has a length of 1, call SetupUDP to move elbow.
if length(prop) == 1
    SetupUDP(SimulationInfo, prop)
    return
end % END IF

% Counts the number of features to be built.
featureCount = 2 * sum(~isnan(prop));

% If there is a value for the last position in prop, subtract 1 from the
% total count. (Only one Thumb joint)
if ~isnan(prop(5))
    featureCount = featureCount - 1;
end

% Initializes the packet
udpPacket = [];
udpPacket = [ udpPacket; int8(featureCount) ];

% Builds the packet with all of the fingers.
udpPacket = IndexUDP(udpPacket, SimulationInfo, prop(1));
udpPacket = MiddleUDP(udpPacket, SimulationInfo, prop(2));
udpPacket = RingUDP(udpPacket, SimulationInfo, prop(3));
udpPacket = LittleUDP(udpPacket, SimulationInfo, prop(4));
udpPacket = ThumbUDP(udpPacket, SimulationInfo, prop(5));


% Sends the packet to MSMS
judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )

end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function udpPacket = ThumbUDP( udpPacket, SimulationInfo, prop )

if isnan(prop)
    return
end

% Parses the ThumbPIP information
compType   = SimulationInfo.Components.ThumbDIP.componentType.charID;
compNumber = SimulationInfo.Components.ThumbDIP.componentNumber.numID;
nDOF       = SimulationInfo.Components.ThumbDIP.nDOF;

% Builds the packet
% [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
udpPacket = [ udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
udpPacket = [ udpPacket; flipud( typecast( int16( compNumber ), 'int8' )' ) ]; % Joint NumID
udpPacket = [ udpPacket; flipud( typecast( int16( nDOF ), 'int8' )' ) ]; % # DOF
m1 = length( udpPacket ) + 1;

% Sets the [m or rad] chunk to 0.
for n = 1:nDOF
    udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; %#ok<AGROW> % zero radians
end % END FOR

% Spans from open to closed state.
rot = prop / 5 * (0.45*pi);

% Inserts rotation values into the packet and sends the packet.
udpPacket( (m1):(m1+3) ) = flipud( typecast( single( rot ), 'int8' )' );

end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function udpPacket = IndexUDP( udpPacket, SimulationInfo, prop )

if isnan(prop)
    return 
end

% Parses the IndexPIP information
compTypePIP   = SimulationInfo.Components.IndexPIP.componentType.charID;
compNumberPIP = SimulationInfo.Components.IndexPIP.componentNumber.numID;
nDOFPIP       = SimulationInfo.Components.IndexPIP.nDOF;

% Parses the IndexDIP information
compTypeDIP   = SimulationInfo.Components.IndexDIP.componentType.charID;
compNumberDIP = SimulationInfo.Components.IndexDIP.componentNumber.numID;
nDOFDIP       = SimulationInfo.Components.IndexDIP.nDOF;

% Builds the packet
% [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
udpPacket = [ udpPacket; int8( compTypePIP(1) ); int8( compTypePIP(2) ) ]; % Joint CharID
udpPacket = [ udpPacket; flipud( typecast( int16( compNumberPIP ), 'int8' )' ) ]; % Joint NumID
udpPacket = [ udpPacket; flipud( typecast( int16( nDOFPIP ), 'int8' )' ) ]; % # DOF
mPIP = length( udpPacket ) + 1;

% Sets the [m or rad] PIP chunk to 0 radians or meters.
udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero radians

udpPacket = [ udpPacket; int8( compTypeDIP(1) ); int8( compTypeDIP(2) ) ]; % Joint CharID
udpPacket = [ udpPacket; flipud( typecast( int16( compNumberDIP ), 'int8' )' ) ]; % Joint NumID
udpPacket = [ udpPacket; flipud( typecast( int16( nDOFDIP ), 'int8' )' ) ]; % # DOF
mDIP = length( udpPacket ) + 1;

% Sets the [m or rad] DIP chunk to 0 radians or meters.
udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero radians

% Spans from open to closed state.
rPIP = prop / 5 * (0.45*pi);
rDIP = prop / 5 * (0.25*pi);

% Inserts rotation values into the packet and sends the packet.
udpPacket( (mPIP):(mPIP+3) ) = flipud( typecast( single( rPIP ), 'int8' )' );
udpPacket( (mDIP):(mDIP+3) ) = flipud( typecast( single( rDIP ), 'int8' )' );

end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function udpPacket = MiddleUDP( udpPacket, SimulationInfo, prop )

if isnan(prop)
    return 
end

% Parses the MiddlePIP information
compTypePIP   = SimulationInfo.Components.MiddlePIP.componentType.charID;
compNumberPIP = SimulationInfo.Components.MiddlePIP.componentNumber.numID;
nDOFPIP       = SimulationInfo.Components.MiddlePIP.nDOF;

% Parses the MiddleDIP information
compTypeDIP   = SimulationInfo.Components.MiddleDIP.componentType.charID;
compNumberDIP = SimulationInfo.Components.MiddleDIP.componentNumber.numID;
nDOFDIP       = SimulationInfo.Components.MiddleDIP.nDOF;

% Builds the packet
% [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
udpPacket = [ udpPacket; int8( compTypePIP(1) ); int8( compTypePIP(2) ) ]; % Joint CharID
udpPacket = [ udpPacket; flipud( typecast( int16( compNumberPIP ), 'int8' )' ) ]; % Joint NumID
udpPacket = [ udpPacket; flipud( typecast( int16( nDOFPIP ), 'int8' )' ) ]; % # DOF
mPIP = length( udpPacket ) + 1;

% Sets the [m or rad] PIP chunk to 0 radians or meters.
udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero radians

udpPacket = [ udpPacket; int8( compTypeDIP(1) ); int8( compTypeDIP(2) ) ]; % Joint CharID
udpPacket = [ udpPacket; flipud( typecast( int16( compNumberDIP ), 'int8' )' ) ]; % Joint NumID
udpPacket = [ udpPacket; flipud( typecast( int16( nDOFDIP ), 'int8' )' ) ]; % # DOF
mDIP = length( udpPacket ) + 1;

% Sets the [m or rad] DIP chunk to 0 radians or meters.
udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero radians

% Spans from open to closed state.
rPIP = prop / 5 * (0.45*pi);
rDIP = prop / 5 * (0.25*pi);

% Inserts rotation values into the packet and sends the packet.
udpPacket( (mPIP):(mPIP+3) ) = flipud( typecast( single( rPIP ), 'int8' )' );
udpPacket( (mDIP):(mDIP+3) ) = flipud( typecast( single( rDIP ), 'int8' )' );
judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )

end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function udpPacket = RingUDP( udpPacket, SimulationInfo, prop )

if isnan(prop)
    return
end

% Parses the RingPIP information
compTypePIP   = SimulationInfo.Components.RingPIP.componentType.charID;
compNumberPIP = SimulationInfo.Components.RingPIP.componentNumber.numID;
nDOFPIP       = SimulationInfo.Components.RingPIP.nDOF;

% Parses the RingDIP information
compTypeDIP   = SimulationInfo.Components.RingDIP.componentType.charID;
compNumberDIP = SimulationInfo.Components.RingDIP.componentNumber.numID;
nDOFDIP       = SimulationInfo.Components.RingDIP.nDOF;

% Builds the packet
% [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
udpPacket = [ udpPacket; int8( compTypePIP(1) ); int8( compTypePIP(2) ) ]; % Joint CharID
udpPacket = [ udpPacket; flipud( typecast( int16( compNumberPIP ), 'int8' )' ) ]; % Joint NumID
udpPacket = [ udpPacket; flipud( typecast( int16( nDOFPIP ), 'int8' )' ) ]; % # DOF
mPIP = length( udpPacket ) + 1;

% Sets the [m or rad] PIP chunk to 0 radians or meters.
udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero radians

udpPacket = [ udpPacket; int8( compTypeDIP(1) ); int8( compTypeDIP(2) ) ]; % Joint CharID
udpPacket = [ udpPacket; flipud( typecast( int16( compNumberDIP ), 'int8' )' ) ]; % Joint NumID
udpPacket = [ udpPacket; flipud( typecast( int16( nDOFDIP ), 'int8' )' ) ]; % # DOF
mDIP = length( udpPacket ) + 1;

% Sets the [m or rad] DIP chunk to 0 radians or meters.
udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero radians

% Spans from open to closed state.
rPIP = prop / 5 * (0.45*pi);
rDIP = prop / 5 * (0.25*pi);

% Inserts rotation values into the packet and sends the packet.
udpPacket( (mPIP):(mPIP+3) ) = flipud( typecast( single( rPIP ), 'int8' )' );
udpPacket( (mDIP):(mDIP+3) ) = flipud( typecast( single( rDIP ), 'int8' )' );

end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function udpPacket = LittleUDP( udpPacket, SimulationInfo, prop )

if isnan(prop)
    return
end

% Parses the LittlePIP information
compTypePIP   = SimulationInfo.Components.LittlePIP.componentType.charID;
compNumberPIP = SimulationInfo.Components.LittlePIP.componentNumber.numID;
nDOFPIP       = SimulationInfo.Components.LittlePIP.nDOF;

% Parses the LittleDIP information
compTypeDIP   = SimulationInfo.Components.LittleDIP.componentType.charID;
compNumberDIP = SimulationInfo.Components.LittleDIP.componentNumber.numID;
nDOFDIP       = SimulationInfo.Components.LittleDIP.nDOF;

% Builds the packet
% [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
udpPacket = [ udpPacket; int8( compTypePIP(1) ); int8( compTypePIP(2) ) ]; % Joint CharID
udpPacket = [ udpPacket; flipud( typecast( int16( compNumberPIP ), 'int8' )' ) ]; % Joint NumID
udpPacket = [ udpPacket; flipud( typecast( int16( nDOFPIP ), 'int8' )' ) ]; % # DOF
mPIP = length( udpPacket ) + 1;

% Sets the [m or rad] PIP chunk to 0 radians or meters.
udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero radians

udpPacket = [ udpPacket; int8( compTypeDIP(1) ); int8( compTypeDIP(2) ) ]; % Joint CharID
udpPacket = [ udpPacket; flipud( typecast( int16( compNumberDIP ), 'int8' )' ) ]; % Joint NumID
udpPacket = [ udpPacket; flipud( typecast( int16( nDOFDIP ), 'int8' )' ) ]; % # DOF
mDIP = length( udpPacket ) + 1;

% Sets the [m or rad] DIP chunk to 0 radians or meters.
udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero radians

% Spans from open to closed state.
rPIP = prop / 5 * (0.45*pi);
rDIP = prop / 5 * (0.25*pi);

% Inserts rotation values into the packet and sends the packet.
udpPacket( (mPIP):(mPIP+3) ) = flipud( typecast( single( rPIP ), 'int8' )' );
udpPacket( (mDIP):(mDIP+3) ) = flipud( typecast( single( rDIP ), 'int8' )' );
end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function SetupUDP( SimulationInfo, prop )

% Parses the ElbowPitch information
jointName  = SimulationInfo.Components.ElbowPitch.Name;
compType   = SimulationInfo.Components.ElbowPitch.componentType.charID;
compNumber = SimulationInfo.Components.ElbowPitch.componentNumber.numID;
seqNum     = SimulationInfo.Components.ElbowPitch.SeqNum;
nDOF       = SimulationInfo.Components.ElbowPitch.nDOF;

% Builds the packet
% [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
udpPacketSize = 1+1*(2+2+2+nDOF*4);
udpPacket = [];
udpPacket = [ udpPacket; int8(1) ]; % change 1 feature
udpPacket = [ udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
udpPacket = [ udpPacket; flipud( typecast( int16( compNumber ), 'int8' )' ) ]; % Joint NumID
udpPacket = [ udpPacket; flipud( typecast( int16( nDOF ), 'int8' )' ) ]; % # DOF
m1 = length( udpPacket ) + 1;

% Sets the [m or rad] chunk to 0.
for n = 1:nDOF
    udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; %#ok<AGROW> % zero radians
end % END FOR

% Warns the user if the packet was built incorrectly. 
if( length( udpPacket ) ~= udpPacketSize )
    warning( 'upd packet size error' ); %#ok<WNTAG>
end % END IF

n = 0.25*(2*pi); % 90 deg in radians

for m = 1:nDOF
    udpPacket( (m1+4*(m-1)):(m1+4*(m-1)+3) ) = flipud( typecast( single( n ), 'int8' )' );
end % END FOR
% Sends the packet to MSMS
judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )

end % END FUNCTION

% EOF