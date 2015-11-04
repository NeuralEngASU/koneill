%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:    MSMSFingerState
% Desc:     Sends the signals to MSMS to control finger joints
% Author:   Kevin O'Neill (Greger Lab)
% Date:     July 10, 2012
% 
% Use:      MSMSFingerState( SimulationInfo, finger, state ) - Will move fingers to an open or closed state
%                           state can be any number. 0 = open
%                                                    != 0 = closed
%           MSMSFingerState( SimulationInfo, 'Index', 1 ) - Closes the Index Proximal-Interphilagial and Index Distal-Interphilagial Joint (IndexPIP/DIP)
%           MSMSFingerState( SimulationInfo, 'Setup', 1 ) - Bends the ElbowPitch joint by 90 deg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
function MSMSFingerState( SimulationInfo, finger, state )
% Moves fingers from an open to a closed state
switch finger
    case 'Thumb'
        ThumbUDP( SimulationInfo, state );
    case 'Index'
        IndexUDP( SimulationInfo, state );
    case 'Middle'
        MiddleUDP( SimulationInfo, state );
    case 'Ring'
        RingUDP( SimulationInfo, state );
    case 'Little'
        LittleUDP( SimulationInfo, state );
    case 'ThumbIndex'
        ThumbUDP( SimulationInfo, state );
        IndexUDP( SimulationInfo, state );
    case 'ThumbMiddle'
        ThumbUDP( SimulationInfo, state );
        MiddleUDP( SimulationInfo, state );
    case 'IndexMiddle'
        IndexUDP( SimulationInfo, state );
        MiddleUDP( SimulationInfo, state );
    case 'TIM'
        ThumbUDP( SimulationInfo, state );
        IndexUDP( SimulationInfo, state );
        MiddleUDP( SimulationInfo, state );
    case 'All'
        ThumbUDP( SimulationInfo, state );
        IndexUDP( SimulationInfo, state );
        MiddleUDP( SimulationInfo, state );
        RingUDP( SimulationInfo, state );
        LittleUDP( SimulationInfo, state );
    case 'Setup'
        SetupUDP( SimulationInfo, state )
    otherwise
end % END SWITCH

end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function ThumbUDP( SimulationInfo, state )
% Parses the ThumbPIP information
jointName  = SimulationInfo.Components.ThumbDIP.Name;
compType   = SimulationInfo.Components.ThumbDIP.componentType.charID;
compNumber = SimulationInfo.Components.ThumbDIP.componentNumber.numID;
seqNum     = SimulationInfo.Components.ThumbDIP.SeqNum;
nDOF       = SimulationInfo.Components.ThumbDIP.nDOF;

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

% Spans from open to closed state.
if logical(state) == true
    rot = (0.45*pi);
else
    rot = 0.0;
end % END IF

% Inserts rotation values into the packet and sends the packet.
udpPacket( (m1):(m1+3) ) = flipud( typecast( single( rot ), 'int8' )' );
judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )

end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function IndexUDP( SimulationInfo, state )
% Parses the IndexPIP information
jointNamePIP  = SimulationInfo.Components.IndexPIP.Name;
compTypePIP   = SimulationInfo.Components.IndexPIP.componentType.charID;
compNumberPIP = SimulationInfo.Components.IndexPIP.componentNumber.numID;
seqNumPIP     = SimulationInfo.Components.IndexPIP.SeqNum;
nDOFPIP       = SimulationInfo.Components.IndexPIP.nDOF;

% Parses the IndexDIP information
jointNameDIP  = SimulationInfo.Components.IndexDIP.Name;
compTypeDIP   = SimulationInfo.Components.IndexDIP.componentType.charID;
compNumberDIP = SimulationInfo.Components.IndexDIP.componentNumber.numID;
seqNumDIP     = SimulationInfo.Components.IndexDIP.SeqNum;
nDOFDIP       = SimulationInfo.Components.IndexDIP.nDOF;

% Builds the packet
% [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
udpPacketSize = 1+1*(2+2+2+nDOFPIP*2*4);
udpPacket = [];
udpPacket = [ udpPacket; int8(2) ]; % change 2 features
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
if logical(state) == true
    rPIP = (0.45*pi);
    rDIP = (0.25*pi);
else
    rPIP = 0.0;
    rDIP = 0.0;
end % END IF

% Inserts rotation values into the packet and sends the packet.
udpPacket( (mPIP):(mPIP+3) ) = flipud( typecast( single( rPIP ), 'int8' )' );
udpPacket( (mDIP):(mDIP+3) ) = flipud( typecast( single( rDIP ), 'int8' )' );
judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )

end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function MiddleUDP( SimulationInfo, state )
% Parses the MiddlePIP information
jointNamePIP  = SimulationInfo.Components.MiddlePIP.Name;
compTypePIP   = SimulationInfo.Components.MiddlePIP.componentType.charID;
compNumberPIP = SimulationInfo.Components.MiddlePIP.componentNumber.numID;
seqNumPIP     = SimulationInfo.Components.MiddlePIP.SeqNum;
nDOFPIP       = SimulationInfo.Components.MiddlePIP.nDOF;

% Parses the MiddleDIP information
jointNameDIP  = SimulationInfo.Components.MiddleDIP.Name;
compTypeDIP   = SimulationInfo.Components.MiddleDIP.componentType.charID;
compNumberDIP = SimulationInfo.Components.MiddleDIP.componentNumber.numID;
seqNumDIP     = SimulationInfo.Components.MiddleDIP.SeqNum;
nDOFDIP       = SimulationInfo.Components.MiddleDIP.nDOF;

% Builds the packet
% [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
udpPacketSize = 1+1*(2+2+2+nDOFPIP*2*4);
udpPacket = [];
udpPacket = [ udpPacket; int8(2) ]; % change 2 features
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
if logical(state) == true
    rPIP = (0.45*pi);
    rDIP = (0.25*pi);
else
    rPIP = 0.0;
    rDIP = 0.0;
end % END IF

% Inserts rotation values into the packet and sends the packet.
udpPacket( (mPIP):(mPIP+3) ) = flipud( typecast( single( rPIP ), 'int8' )' );
udpPacket( (mDIP):(mDIP+3) ) = flipud( typecast( single( rDIP ), 'int8' )' );
judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )

end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function RingUDP( SimulationInfo, state )
% Parses the RingPIP information
jointNamePIP  = SimulationInfo.Components.RingPIP.Name;
compTypePIP   = SimulationInfo.Components.RingPIP.componentType.charID;
compNumberPIP = SimulationInfo.Components.RingPIP.componentNumber.numID;
seqNumPIP     = SimulationInfo.Components.RingPIP.SeqNum;
nDOFPIP       = SimulationInfo.Components.RingPIP.nDOF;

% Parses the RingDIP information
jointNameDIP  = SimulationInfo.Components.RingDIP.Name;
compTypeDIP   = SimulationInfo.Components.RingDIP.componentType.charID;
compNumberDIP = SimulationInfo.Components.RingDIP.componentNumber.numID;
seqNumDIP     = SimulationInfo.Components.RingDIP.SeqNum;
nDOFDIP       = SimulationInfo.Components.RingDIP.nDOF;

% Builds the packet
% [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
udpPacketSize = 1+1*(2+2+2+nDOFPIP*2*4);
udpPacket = [];
udpPacket = [ udpPacket; int8(2) ]; % change 2 features
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
if logical(state) == true
    rPIP = (0.45*pi);
    rDIP = (0.25*pi);
else
    rPIP = 0.0;
    rDIP = 0.0;
end % END IF

% Inserts rotation values into the packet and sends the packet.
udpPacket( (mPIP):(mPIP+3) ) = flipud( typecast( single( rPIP ), 'int8' )' );
udpPacket( (mDIP):(mDIP+3) ) = flipud( typecast( single( rDIP ), 'int8' )' );
judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )

end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function LittleUDP( SimulationInfo, state )
% Parses the LittlePIP information
jointNamePIP  = SimulationInfo.Components.LittlePIP.Name;
compTypePIP   = SimulationInfo.Components.LittlePIP.componentType.charID;
compNumberPIP = SimulationInfo.Components.LittlePIP.componentNumber.numID;
seqNumPIP     = SimulationInfo.Components.LittlePIP.SeqNum;
nDOFPIP       = SimulationInfo.Components.LittlePIP.nDOF;

% Parses the LittleDIP information
jointNameDIP  = SimulationInfo.Components.LittleDIP.Name;
compTypeDIP   = SimulationInfo.Components.LittleDIP.componentType.charID;
compNumberDIP = SimulationInfo.Components.LittleDIP.componentNumber.numID;
seqNumDIP     = SimulationInfo.Components.LittleDIP.SeqNum;
nDOFDIP       = SimulationInfo.Components.LittleDIP.nDOF;

% Builds the packet
% [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
udpPacketSize = 1+1*(2+2+2+nDOFPIP*2*4);
udpPacket = [];
udpPacket = [ udpPacket; int8(2) ]; % change 2 features
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
if logical(state) == true
    rPIP = (0.45*pi);
    rDIP = (0.25*pi);
else
    rPIP = 0.0;
    rDIP = 0.0;
end % END IF

% Inserts rotation values into the packet and sends the packet.
udpPacket( (mPIP):(mPIP+3) ) = flipud( typecast( single( rPIP ), 'int8' )' );
udpPacket( (mDIP):(mDIP+3) ) = flipud( typecast( single( rDIP ), 'int8' )' );
judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )
end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function SetupUDP( SimulationInfo, state )

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