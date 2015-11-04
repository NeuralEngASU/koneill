% Handle Vicon data stream

%% Setup

dllLoc = '';
NET.addAssembly(dllLoc);

vicon = NetThread.init();
Vicon.SetIP('128.0.0.1');
Vicon.SetPort(9090);

Vicon.StartThread();

points = zeros(26,3); % 3 coords for each marker.
angles = zeros(7,4);  % Thumb, Index, Middle, Ring, Little, Right Orient, Left Orient
loc = zeros(2,3);     % Right Loc, Left Loc

for ii = 1:100
%% Get Points
points = Vicon.GetPoints();

%% Process Data
[angles, loc] = Points2Angles(points);

%% Send Angles

end %END FOR

% EOF