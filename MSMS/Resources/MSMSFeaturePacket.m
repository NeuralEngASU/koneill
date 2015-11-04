jointName = SimulationInfo.Components.ElbowPitch.Name;
componentType = SimulationInfo.Components.ElbowPitch.componentType.charID;
componentNumber = SimulationInfo.Components.ElbowPitch.componentNumber.numID;
seqNum = SimulationInfo.Components.ElbowPitch.SeqNum;
nDOF = 1;

% jointName = 'PART Palm assy_1__INDEX_METACAR_Ver 3_1';
% componentType = ' J';
% componentNumber = 1148;
% seqNum = 18;
% nDOF = 1;
% 
% jointName = 'SphereJoint';
% componentType = ' J';
% componentNumber = 200;
% seqNum = 37;
% nDOF = 3;
% 
% jointName = 'SternoClavicular';
% componentType = ' J';
% componentNumber = 39;
% seqNum = 4;
% nDOF = 2;

% jointName = 'CylinderJoint';
% componentType = ' J';
% componentNumber = 207;
% seqNum = 38;
% nDOF = 3;

% jointName = 'PART Thumb proximal assy_1__PART Thumb distal assy_1';
% componentType = ' J';
% componentNumber = 1147;
% seqNum = 36;
% nDOF = 1;

% jointName = 'PART Middle proximal_1__PART Middle mid assy_1';
% componentType = ' J';
% componentNumber = 1149;
% seqNum = 23;
% nDOF = 1;

% jointName = 'ChairJoint';
% componentType = ' J';
% componentNumber = 300;
% seqNum = -1;
% nDOF = 3;

% jointName = 'GlenoHumeral';
% componentType = ' J';
% componentNumber = 36;
% seqNum = 5;
% nDOF = 3;

% jointName = 'Neck';
% componentType = ' J';
% componentNumber = 41;
% seqNum = 3;
% nDOF = 3;

% jointName = 'INDEX_METACAR_Ver 3_1__PART Index proximal assy_1';
% componentType = ' J';
% componentNumber = 1169;
% seqNum = 19;
% nDOF = 1;

% jointName = 'THUMB_YOKE_4DEG Version 2_1__PART thumb base assy_1';
% componentType = ' J';
% componentNumber = 1151;
% seqNum = 34;
% nDOF = 1;

% jointName = 'SHOULDER GEARMOTOR ASSEMBLY_1__SHOULDER TO HR_1';
% componentType = ' J';
% componentNumber = 1158;
% seqNum = 12;
% nDOF = 1;

% jointName = 'PART little proximal assy_1__PART little mid assy_1';
% componentType = ' J';
% componentNumber = 1155;
% seqNum = 31;
% nDOF = 1;

% jointName = 'COBOT MAIN ASSEMBLY_1__WRIST_1';
% componentType = ' J';
% componentNumber = 1146;
% seqNum = 16;
% nDOF = 1;

% jointName = 'WRIST_1__PART Palm assy_1';
% componentType = ' J';
% componentNumber = 1153;
% seqNum = 17;
% nDOF = 1;

% jointName = 'SphereJoint';
% componentType = ' J';
% componentNumber = 200;
% seqNum = 37;
% nDOF = 3;

% jointName = 'CubeJoint';
% componentType = ' J';
% componentNumber = 301;
% seqNum = 39;
% nDOF = 3;

% jointName = 'TableJoint';
% componentType = ' J';
% componentNumber = 92;
% seqNum = -1;
% nDOF = 3;

% jointName = 'SternoClavicularL';
% componentType = ' J';
% componentNumber = 40;
% seqNum = 9;
% nDOF = 2;

% sending packet to MSMS
udpPacketSize = 1+1*(2+2+2+nDOF*4);
udpPacket = [];
udpPacket = [ udpPacket; int8(1) ]; % change 1 feature
udpPacket = [ udpPacket; int8( componentType(1) ); int8( componentType(2) ) ]; % Joint ID
udpPacket = [ udpPacket; flipud( typecast( int16( componentNumber ), 'int8' )' ) ]; % Neck ID
udpPacket = [ udpPacket; flipud( typecast( int16( nDOF ), 'int8' )' ) ]; % 3 DOF in radians
m1 = length( udpPacket ) + 1;
for n = 1:nDOF
    udpPacket = [ udpPacket; flipud( typecast( single( 0 ), 'int8' )' ) ]; %#ok<AGROW> % zero radians
end
if( length( udpPacket ) ~= udpPacketSize )
    warning( 'upd packet size error' ); %#ok<WNTAG>
end

for n = 0.0:(0.025*(2*pi)):(2*pi)
    for m = 1:nDOF
        udpPacket( (m1+4*(m-1)):(m1+4*(m-1)+3) ) = flipud( typecast( single( n ), 'int8' )' );
    end
    judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )
    pause( 0.01 );
end
fprintf( '\n' );

