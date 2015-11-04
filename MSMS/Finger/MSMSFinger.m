%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:    MSMSFinger
% Desc:     Sends the signals to MSMS to control finger joints
% Author:   Kevin O'Neill (Greger Lab)
% Date:     July 2, 2012
% 
% Use:      MSMSFinger( SimulationInfo, finger, state ) - Will move fingers to an open or closed state
%           MSMSFinger( SimulationInfo, 'Index', 1 ) - Closes the Index Proximal-Interphilagial Joint (IndexPIP)
%           MSMSFinger( SimulationInfo, 'Setup', 1 ) - Bends the ElbowPitch joint by 90 deg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
function MSMSFinger( SimulationInfo, finger, state )
% Moves fingers from an open to a closed state

% SimulationInfo = ParseXML('D:\CodeRepo\PatientCart\MATLAB\MSMS\simulationSetup.xml');

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
function IndexUDP( SimulationInfo, state )
% Parses the IndexPIP information
jointName  = SimulationInfo.Components.IndexPIP.Name;
compType   = SimulationInfo.Components.IndexPIP.componentType.charID;
compNumber = SimulationInfo.Components.IndexPIP.componentNumber.numID;
seqNum     = SimulationInfo.Components.IndexPIP.SeqNum;
nDOF       = SimulationInfo.Components.IndexPIP.nDOF;

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
    for n = 0.0:(0.025*(2*pi)):(.25*pi)
        for m = 1:nDOF
            udpPacket( (m1+4*(m-1)):(m1+4*(m-1)+3) ) = flipud( typecast( single( n ), 'int8' )' );
        end % END FOR
        judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )
        pause( 0.01 );
    end % END FOR
else
    for n = (0.25*pi):-(0.025*(2*pi)):0.0
        for m = 1:nDOF
            udpPacket( (m1+4*(m-1)):(m1+4*(m-1)+3) ) = flipud( typecast( single( n ), 'int8' )' );
        end
        judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )
        pause( 0.01 );
    end
end % END IF

end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function MiddleUDP( SimulationInfo, state )
% Parses the MiddlePIP information
jointName  = SimulationInfo.Components.MiddlePIP.Name;
compType   = SimulationInfo.Components.MiddlePIP.componentType.charID;
compNumber = SimulationInfo.Components.MiddlePIP.componentNumber.numID;
seqNum     = SimulationInfo.Components.MiddlePIP.SeqNum;
nDOF       = SimulationInfo.Components.MiddlePIP.nDOF;

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
    for n = 0.0:(0.025*(2*pi)):(.25*pi)
        for m = 1:nDOF
            udpPacket( (m1+4*(m-1)):(m1+4*(m-1)+3) ) = flipud( typecast( single( n ), 'int8' )' );
        end % END FOR
        judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )
        pause( 0.01 );
    end % END FOR
else
    for n = (0.25*pi):-(0.025*(2*pi)):0.0
        for m = 1:nDOF
            udpPacket( (m1+4*(m-1)):(m1+4*(m-1)+3) ) = flipud( typecast( single( n ), 'int8' )' );
        end % END FOR
        judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )
        pause( 0.01 );
    end % END FOR
end % END IF

end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function RingUDP( SimulationInfo, state )
% Parses the RingPIP information
jointName  = SimulationInfo.Components.RingPIP.Name;
compType   = SimulationInfo.Components.RingPIP.componentType.charID;
compNumber = SimulationInfo.Components.RingPIP.componentNumber.numID;
seqNum     = SimulationInfo.Components.RingPIP.SeqNum;
nDOF       = SimulationInfo.Components.RingPIP.nDOF;

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
    for n = 0.0:(0.025*(2*pi)):(.25*pi)
        for m = 1:nDOF
            udpPacket( (m1+4*(m-1)):(m1+4*(m-1)+3) ) = flipud( typecast( single( n ), 'int8' )' );
        end % END FOR
        judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )
        pause( 0.01 );
    end % END FOR
else
    for n = (0.25*pi):-(0.025*(2*pi)):0.0
        for m = 1:nDOF
            udpPacket( (m1+4*(m-1)):(m1+4*(m-1)+3) ) = flipud( typecast( single( n ), 'int8' )' );
        end % END FOR 
        judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) ) % Sends packet to MSMS
        pause( 0.01 );
    end % END FOR
end % END IF

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
    for n = 0.0:(0.025*(2*pi)):(.25*pi)
        for m = 1:nDOF
            udpPacket( (m1+4*(m-1)):(m1+4*(m-1)+3) ) = flipud( typecast( single( n ), 'int8' )' );
        end % END FOR
        judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )
        pause( 0.01 );
    end % END FOR
else
    for n = (0.25*pi):-(0.025*(2*pi)):0.0
        for m = 1:nDOF
            udpPacket( (m1+4*(m-1)):(m1+4*(m-1)+3) ) = flipud( typecast( single( n ), 'int8' )' );
        end % END FOR 
        judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) ) % Sends packet to MSMS
        pause( 0.01 );
    end % END FOR
end % END IF

end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function LittleUDP( SimulationInfo, state )
% Parses the LittlePIP information
jointName  = SimulationInfo.Components.LittlePIP.Name;
compType   = SimulationInfo.Components.LittlePIP.componentType.charID;
compNumber = SimulationInfo.Components.LittlePIP.componentNumber.numID;
seqNum     = SimulationInfo.Components.LittlePIP.SeqNum;
nDOF       = SimulationInfo.Components.LittlePIP.nDOF;

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
    for n = 0.0:(0.025*(2*pi)):(.25*pi)
        for m = 1:nDOF
            udpPacket( (m1+4*(m-1)):(m1+4*(m-1)+3) ) = flipud( typecast( single( n ), 'int8' )' );
        end % END FOR
        judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )
        pause( 0.01 );
    end % END FOR
else
    for n = (0.25*pi):-(0.025*(2*pi)):0.0
        for m = 1:nDOF
            udpPacket( (m1+4*(m-1)):(m1+4*(m-1)+3) ) = flipud( typecast( single( n ), 'int8' )' );
        end % END FOR 
        judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) ) % Sends packet to MSMS
        pause( 0.01 );
    end % END FOR
end % END IF

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