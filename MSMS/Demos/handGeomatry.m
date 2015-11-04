% relPos = [0.4503;-0.1863;-0.0258];
% 
% % Palmpos = 0.4503, -0.1863, -0.0258 m
% % Palm vector = 0,0,-1
% % Palm rot = 0, 330, 90 deg
% % http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/
% 
% for n = 1 : 4
%     relPos(:,n+1) = relPos(:,n) + PlacementMat(1:3,n);
% end
% 
%  relPos(4,:) = [ (-30/360)*(2*pi), 0, (10/360)*(2*pi), (70/360)*(2*pi), (30/360)*(2*pi)];
% relpos2 = [];

% relPos              = zeros(7,6); % Initialize matrix
% relPos(1:3,2:end-1) = fliplr(PlacementMat(1:3,:)); % Order matrix from distal to proximal (with the desired position at the first column)
% relPos(4,:)         = [ 0, (30/360)*(2*pi), (70/360)*(2*pi), (10/360)*(2*pi), 0, (-30/360)*(2*pi)]; % Degrees of the joints from distal to proximal
% relPos(1:3,1)         = [-0.0095, -0.0140, 0.0055]; % Offset of finger 'pad' from distal joint (guess: 0.5 inch = 0.0127 m)
% relPos(5:7,2:end-1) = fliplr(PlacementMat(4:6,:));
% relPos(6,6)         = -1;

%%

% load('D:\CodeRepo\PatientCart\MATLAB\MSMS\Demos\relPos.mat');

% load('C:\CodeRepo\Lab\PatientCart\MATLAB\MSMS\Demos\relPosIndex.mat');

% load('C:\CodeRepo\Lab\PatientCart\MATLAB\MSMS\Demos\relPosIndexOld0718.mat')

load('C:\CodeRepo\Lab\PatientCart\MSMS\PositionMat\relPos.mat')

relPos = relPos{2};

% relPos = relPosIndex;

% load('C:\CodeRepo\Lab\PatientCart\MATLAB\MSMS\Demos\relPosMiddle.mat')
% Removed the palm offset.
% load('C:\CodeRepo\Lab\PatientCart\MATLAB\MSMS\DemosOld\TESTrelPos.mat')

% relPos = fliplr(relPosThumb);
% 
relPos = fliplr(relPos);
for n = 2: size(relPos,2) %8
    relPos(1:3,n) = relPos(1:3,n) + relPos(1:3,n-1);
end
relPos = fliplr(relPos);
% 
% relPos = fliplr(relPos);
% 
% relPosThumb = relPos
% 
% relPos(:,end+1) = relPos(:,end);
% relPos(:,end-1) = [0;0;0;pi/2;0;0;1];
% 
% offsetNew = relPos(1:3, 1);

%%

% % Index/Middle
rot3(1,:) = linspace(0,(10/360)*(2*pi),100);
rot3(2,:) = linspace(0,(70/360)*(2*pi),100);
rot3(3,:) = linspace(0,(30/360)*(2*pi),100);

% Thumb
% rot3(1,:) = linspace(0,(45/360)*(2*pi),100);
% rot3(2,:) = linspace(0,(10/360)*(2*pi),100);

pos3d = [0;0;0];

for j = 1:100
    
%     relPos(4,2:3) = [ rot3(2,j), rot3(1,j)];%, rot3(2,j), rot3(1,j)];
    relPos(4,2:4) = [ rot3(3,j), rot3(2,j), rot3(1,j) ];

    
    offsetNew = relPos(1:3, 1);
    
    for k = 1:size(relPos,2)-2%7
        
        theta = relPos(  4, k+1);
        abc   = relPos(1:3, k+1);
        xyz   = offsetNew;
        uvw   = relPos(5:7, k+1);
        
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
        
        
        %     offsetNew = xyzNew + abc;
        offsetNew= xyzNew;
        
    end
    pos3d(:,j) = xyzNew;
    
    plot3(pos3d(1,:),pos3d(2,:),pos3d(3,:))

%     plot(pos3d(2,:),pos3d(3,:))
%     axis([-0.5,00.1,0,-.1,-0.2,0,0.035])
    axis equal
    drawnow
%     xlabel('y')
%     ylabel('z')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    pause(0.01)
end

hold on
plot3(pos3d(1,1),pos3d(2,1),pos3d(3,1), 'og')
plot3(pos3d(1,end),pos3d(2,end),pos3d(3,end), 'or')
hold off

%%
% 
% absPosEnd(1) = relPos(1,end) + pos3d(1,end);
% absPosEnd(2) = relPos(2,end) + pos3d(2,end);
% absPosEnd(3) = relPos(3,end) + pos3d(3,end);

% absPosEnd = relPos(1:3,end) + pos3d(:,end);
absPosEnd = pos3d(:,end);

% position = [ 0.5750; -0.1435; -0.0110 ];
position = [ 0.5850; -0.1610; 0.0 ];

difference = (position - absPosEnd);

sqrt(sum(difference.^2))


% palmPos = [0.4503, -0.1863, -0.0258];
% palmPos(1) = relPos(1, end)-offsetNew(2);
% palmPos(2) = relPos(2, end)+offsetNew(1);
% palmPos(3) = relPos(3, end)+offsetNew(3);
%       

% for k=2%1:4
%     
%     theta = relPosCopy(4,k);
%     abc   = relPos(1:3,k);
%     xyz   = [relPos(1:3,k+1);1];
%     u     = relPosCopy(5,k);
%     v     = relPosCopy(6,k);
%     w     = relPosCopy(7,k);
%     
%     Tp = eye(4);
%     Tp(1:3,4) = abc;
%     Rz = [cos(theta), -sin(theta), 0, 0; ...
%         sin(theta),  cos(theta), 0, 0; ...
%         0,           0, 1, 0; ...
%         0,           0, 0, 1];
%     Txz = [u/sqrt(u^2 + v^2), v/sqrt(u^2 + v^2), 0, 0; ...
%         -v/sqrt(u^2 + v^2), u/sqrt(u^2 + v^2), 0, 0; ...
%         0,                 0, 1, 0; ...
%         0,                 0, 0, 1];
%     Tz = [         w/sqrt(u^2+v^2+w^2), 0, -sqrt(u^2+v^2)/sqrt(u^2+v^2+w^2), 0; ...
%                                         0, 1,                             0, 0; ...
%           sqrt(u^2+v^2)/sqrt(u^2+v^2+w^2), 0,           w/sqrt(u^2+v^2+w^2), 0; ...
%                                         0, 0,                             0, 1];
%     
%     
%     
%     
%     
%     tmp = (inv(Tp) * inv(Txz) * inv(Tz) * Rz * Tz * Txz * Tp) * xyz;
%     relPos(1:3,k+1) = tmp(1:3)
% 
% end

% relPos = relPosCopy;

% Final Point + some offset (offset being the point of the pad on the
% finger.)



