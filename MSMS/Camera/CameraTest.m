for compNumber = 118%90:120
    
    compType   = 'HT';
    nDOF = 6;
    
    % Builds the packet
    % [# Features Changed][Component CharID][Component NumID][# DOF][ m or rad]
    udpPacket = [];
    udpPacket = [ udpPacket; int8(1) ]; % change 1 feature
    udpPacket = [ udpPacket; int8( compType(1) ); int8( compType(2) ) ]; % Joint CharID
    udpPacket = [ udpPacket; flipud( typecast( int16( compNumber ), 'int8' )' ) ]; % Joint NumID
    udpPacket = [ udpPacket; flipud( typecast( int16( 14 ), 'int8' )' ) ]; % FEATURE TO CHANGE 14 = Head Tracking aka Camera
    m1 = length( udpPacket ) + 1;
    
    % Sets the [m or rad] chunk to 0.
    for n = 1:nDOF
        udpPacket = [ udpPacket; flipud( typecast( single( 2 ), 'int8' )' ) ]; %#ok<AGROW> % zero radians
    end % END FOR
    
    % % Spans from open to closed state.
    % if logical(state) == true
    %     for n = 0.0:(0.025*(2*pi)):(.25*pi)
    %         for m = 1:nDOF
    %             udpPacket( (m1+4*(m-1)):(m1+4*(m-1)+3) ) = flipud( typecast( single( n ), 'int8' )' );
    %         end % END FOR
    %         judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )
    %         pause( 0.01 );
    %     end % END FOR
    % else
    %     for n = (0.25*pi):-(0.025*(2*pi)):0.0
    %         for m = 1:nDOF
    %             udpPacket( (m1+4*(m-1)):(m1+4*(m-1)+3) ) = flipud( typecast( single( n ), 'int8' )' );
    %         end
    %         judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )
    %         pause( 0.01 );
    %     end
    % end % END IF
    
    judp( 'send', 11114, '127.0.0.1', typecast( udpPacket, 'int8' ) )
    
    compNumber
    pause(0.05)
    
    
end
