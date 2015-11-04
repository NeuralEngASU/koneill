%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:    ObjectControl
% Desc:     Sends information to MSMS via UDP to control on screen objects.
% Author:   Kevin O'Neill (Greger Lab)
% Date:     July 11, 2012
% 
% Use:      
%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function ObjectControl( SimInfo, prop )

% Counts the number of features to be built.
featureCount = sum(logical(sum(~isnan(prop()))));

% Initializes the packet
udpPacket = [];
udpPacket = [ udpPacket; int8(3) ];

% Builds the packet with all of the fingers.
udpPacket = CubeUDP(udpPacket, SimInfo, prop(:,1));
udpPacket = CylinderUDP(udpPacket, SimInfo, prop(:,2));
udpPacket = SphereUDP(udpPacket, SimInfo, prop(:,3));

% Sends the packet to MSMS
judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )
% judp( 'send', 11114, '155.101.184.12', typecast( udpPacket, 'int8' ) )

end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function udpPacket = CubeUDP( udpPacket, SimInfo, prop )

if ~sum(~isnan(prop))
    return 
end

% Parses the IndexPIP information
compType   = SimInfo.Components.Cube.componentType.charID;
compNumber = SimInfo.Components.Cube.componentNumber.numID;
nDOF       = SimInfo.Components.Cube.nDOF;

% Builds the packet
% [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
udpPacket = [ udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
udpPacket = [ udpPacket; flipud( typecast( int16( compNumber ), 'int8' )' ) ]; % Joint NumID
udpPacket = [ udpPacket; flipud( typecast( int16( nDOF ), 'int8' )' ) ]; % # DOF
messLen = length( udpPacket ) + 1;

% Sets the [m or rad] PIP chunk to 0 radians or meters.
for k = 1 : nDOF
    udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; % Zero radians
end % END FOR

% Inserts rotation values into the packet and sends the packet.
for m = 1 : nDOF
    udpPacket( (messLen+4*(m-1)):(messLen+4*(m-1)+3) ) = flipud( typecast( single( prop(m) ), 'int8' )' );
end % END FOR

end % END FUNCTION

% EOF