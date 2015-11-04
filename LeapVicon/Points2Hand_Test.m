%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Points2Hand
%   Author: Kevin O'Neill
%   Date: 2015/02/27
%   Desc: 
%       This file contains pseudo-code for how to turn points into hand
%       angles. There are 27 DOFs that must be calculated to recreate hand
%       motion. Reduced and constrained models can have 23 DOFs to
%       calculate (excludes DIP since DIP and IIP are usually the same).
%
%       Finger Joint Angles (16):
%           Adb (Abduction / Aduction)
%           PIP (Proximal Interphalangeal joint)
%           IIP (Intermediate Interphalangeal joint)
%           DIP (Distal Interphalangeal joint)
%
%       Thumb Joint Angles (4):
%           Adb (Abduction / Aduction)
%           CMC (Carpometacarpal joint)            
%           MP (Metacarpophalangeal  joint)
%           IP (Interphalangeal joint)
%
%       Wrist Joint Angles(1):
%           WF (Wrist flexion/extension)sleep
%
%       Palm Degrees of Freedom (6):
%           Translation (3)
%           Rotation (3)
%
%       Marker Indicies: 23
%           One marker on each finger nail: 5
%           One marker on each finger joint: 15
%           One marker, centered, on the back of the hand: 1
%           One marker on the wrist: 1
%           One marker on the dorsal surface of the forearm: 1
%
%       Labels will follow a simplistic pattern for ease of typing and
%       reading. Each label will be a single set of (x,y,z) points in space.
%
%       Index Knuckle ( I0 )      - PIP and Adb Joint
%       Index IIP ( I1 )          - IIP Joint
%       Index DIP ( I2 )          - DIP Joint
%       Index Finger Nail ( I3 )  - Finger Tip
%
%       Repeat for each finger and thumb (T,I,M,R,L)       
%
%       Hand Back ( HB )
%
%       Wrist Joint ( WJ )
%
%       Forearm ( FA )
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load and seperate data

rotMatx = @(x) [      1,       0,        0;...
                      0, cosd(x), -sind(x);...
                      0, sind(x),  cosd(x)];
                  
rotMaty = @(x) [cosd(x),       0,  sind(x);...
                      0,       1,        0;...
               -sind(x),       0,  cosd(x)];
           
rotMatz = @(x) [cosd(x),-sind(x),        0;...
                sind(x), cosd(x),        0;...
                      0,       0,        1];

jointAngles = [55,   45,   90,  90;
               -10,   45,   80,   90;...
               0,   45,   70,   90;...
               10,   45,   50,   90;...
               20,   45,   45,   90];
    
jointAngles(1,1) = jointAngles(1,1) -75.9638;            

handAngles = [10, 0, 45];

ThumbMarkers = zeros(4,3);
IndexMarkers = zeros(4,3);
MiddleMarkers = zeros(4,3);
RingMarkers = zeros(4,3);
LittleMarkers = zeros(4,3);

translation = [5,-6,7];
handRotMat = rotMatz(handAngles(2)) * rotMaty(handAngles(3)) * rotMatx(handAngles(1));

HB = [0;0;0]' + translation;
WF = HB + [0,-1,0]*handRotMat;

for i = 1:4

    if i == 1
        ThumbMarkers(1,:)  = [-1.5, -0.5 , 0] * handRotMat + HB;
        IndexMarkers(1,:)  = [-1  ,  1   , 0] * handRotMat + HB;
        MiddleMarkers(1,:) = [ 0  ,  1   , 0] * handRotMat + HB;
        RingMarkers(1,:)   = [ 1  ,  0.95, 0] * handRotMat + HB;
        LittleMarkers(1,:) = [ 2  ,  0.80, 0] * handRotMat + HB;
%     elseif i == 2
%         ThumbMarkers(i,:)  = (rotMatz(jointAngles(1,1))*rotMatx(sum(jointAngles(1,2:i)))*handRotMat* [0;-1;0] .* [1;-1;1] + ThumbMarkers(i-1,:)')';
%         IndexMarkers(i,:)  = (rotMatz(jointAngles(2,1))*rotMatx(sum(jointAngles(2,2:i))) *handRotMat * [0;-1;0] .* [1;-1;1]+ IndexMarkers(i-1,:)')';
%         MiddleMarkers(i,:) = (rotMatz(jointAngles(3,1))*rotMatx(sum(jointAngles(3,2:i)))*handRotMat * [0;-1;0] .* [1;-1;1] + MiddleMarkers(i-1,:)')';
%         RingMarkers(i,:)   = (rotMatz(jointAngles(4,1))*rotMatx(sum(jointAngles(4,2:i)))*handRotMat * [0;-1;0] .* [1;-1;1] + RingMarkers(i-1,:)')';
%         LittleMarkers(i,:) = (rotMatz(jointAngles(5,1))*rotMatx(sum(jointAngles(5,2:i)))*handRotMat * [0;-1;0] .* [1;-1;1] + LittleMarkers(i-1,:)')';
       
    else
        ThumbMarkers(i,:)  = (rotMatz(jointAngles(1,1))*rotMatx(sum(jointAngles(1,2:i))) * [0;-1;0] .* [1;-1;1] + ThumbMarkers(i-1,:)')';
        IndexMarkers(i,:)  = (rotMatz(jointAngles(2,1))*rotMatx(sum(jointAngles(2,2:i))) * [0;-1;0] .* [1;-1;1] + IndexMarkers(i-1,:)')';
        MiddleMarkers(i,:) = (rotMatz(jointAngles(3,1))*rotMatx(sum(jointAngles(3,2:i))) * [0;-1;0] .* [1;-1;1] + MiddleMarkers(i-1,:)')';
        RingMarkers(i,:)   = (rotMatz(jointAngles(4,1))*rotMatx(sum(jointAngles(4,2:i))) * [0;-1;0] .* [1;-1;1] + RingMarkers(i-1,:)')';
        LittleMarkers(i,:) = (rotMatz(jointAngles(5,1))*rotMatx(sum(jointAngles(5,2:i))) * [0;-1;0] .* [1;-1;1] + LittleMarkers(i-1,:)')';
        
    end % END IF
    
end % END FOR

plotOn = 1;


if plotOn
    figure(1)
    hold on
    
    scatter3(HB(1),HB(2),HB(3), 'k');
    plot3([HB(1), ThumbMarkers(1,1)], [HB(2), ThumbMarkers(1,2)], [HB(3),ThumbMarkers(1,3)], 'k')
    plot3([HB(1), IndexMarkers(1,1)], [HB(2), IndexMarkers(1,2)], [HB(3),IndexMarkers(1,3)], 'k')
    plot3([HB(1), MiddleMarkers(1,1)], [HB(2), MiddleMarkers(1,2)], [HB(3),MiddleMarkers(1,3)], 'k')
    plot3([HB(1), RingMarkers(1,1)], [HB(2), RingMarkers(1,2)], [HB(3),RingMarkers(1,3)], 'k')
    plot3([HB(1), LittleMarkers(1,1)], [HB(2), LittleMarkers(1,2)], [HB(3),LittleMarkers(1,3)], 'k')
    
    scatter3(ThumbMarkers(:,1), ThumbMarkers(:,2), ThumbMarkers(:,3), 'r');
    plot3(ThumbMarkers(:,1), ThumbMarkers(:,2), ThumbMarkers(:,3), 'r')
    
    scatter3(IndexMarkers(:,1), IndexMarkers(:,2), IndexMarkers(:,3), 'g');
    plot3(IndexMarkers(:,1), IndexMarkers(:,2), IndexMarkers(:,3), 'g')
    
    scatter3(MiddleMarkers(:,1), MiddleMarkers(:,2), MiddleMarkers(:,3), 'b');
    plot3(MiddleMarkers(:,1), MiddleMarkers(:,2), MiddleMarkers(:,3), 'b')
    
    scatter3(RingMarkers(:,1), RingMarkers(:,2), RingMarkers(:,3), 'm');
    plot3(RingMarkers(:,1), RingMarkers(:,2), RingMarkers(:,3), 'm')
    
    scatter3(LittleMarkers(:,1), LittleMarkers(:,2), LittleMarkers(:,3), 'c');
    plot3(LittleMarkers(:,1), LittleMarkers(:,2), LittleMarkers(:,3), 'c')
    hold off

end % END IF
%% Important variables anf functions

% Define a function [vecMag] such that it returns the magnitude of a matrix
% [x] of vectors [x_i].
vecMag = @(x) sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2);

% Rotation matrix decomposition
rotAngles = @(x) [(atan2(x(3,2,:), x(3,3,:)) .* (180/pi)),...
                  -acos(x(3,1,:)).*(180/pi),...  %(atan2(-x(3,1,:), sqrt(x(3,2,:).^2 + x(3,3,:).^2)) .* (180/pi)),...
                  (atan2(x(2,1,:), x(1,1,:)) .* (180/pi))];

              
% Allocate variable to store angles
calcAngles = zeros(5,4);

% Create FingerMarkers matrix
FingerMarkers = [ThumbMarkers;IndexMarkers;MiddleMarkers;RingMarkers;LittleMarkers];

% Calculate the vectors from one joint to another
vecFingers = FingerMarkers(2:end,:) - FingerMarkers(1:end-1,:);
vecFingers(~mod(1:length(vecFingers), 4)',:) = []; % Cull excess vectors (eg: IndexTip to Middle Knuckle)

% Rotation matrix for each joint
rotMat = zeros(3,3,5*3); % One rotation matrix for each joint for each finger. Thumb MC, IIP, PIP, Index MC,...

%% Finger Angles (Not Knuckle)

tic
% Create an index variable 
idx1 = 1:length(vecFingers);
idx1(~mod(idx1,3)) = [];
idx2 = idx1+1;

% Allocate variable to store calculated angles between distal bones of the
% fingers. (2 joints for 5 fingers -> 10 indicies)
anglesTmp = zeros(1,10);

%  Calculate the angle for each distal and intermediate joint of the finger
% theta = acosd( dot(v,w) / (||v|| * ||w||));
anglesTmp = acosd( dot(vecFingers(idx1,:),vecFingers(idx2,:),2) ./ (vecMag(vecFingers(idx1,:)).*vecMag(vecFingers(idx2,:))));

% Store the calculated angles in the [calcAngles] matrix. Use the [real]
% values as imaginary values occur when the angle is 0 (or 180);
calcAngles(:,3:4) = real(reshape(anglesTmp',2,5)');


%% Rotation Matrix decomposition: Tait-Bryan Angles
% Bone basis orientation:
%   x-axis: pointing 'down' out of the bone if the hand is palm down.
%   y-axis: point 'right' out of the bone if the hand is palm down.
%   z-axis: along bone from proximal end to distal end.
              

% FingerMarkers = FingerMarkers + repmat([5,6,7], length(FingerMarkers),1);
% Calculate the vectors from one joint to another


% Create hand orientation basis
% Wrist to knuckles (z-axis)
% WF = HB' + [0,-1,0];
handz = HB-WF;
handz = handz./norm(handz);

% Middle knuckle and Index knuckle, respectivly
handx = cross(HB-FingerMarkers(5,:), HB-FingerMarkers(9,:));
handx = handx./norm(handx);

handy = cross(handx,handz);
handy = handy./norm(handy);

handRot = [handx', handy', handz'];

handRotDecom = rotAngles(handRot);

% Create Normals between bones. Ex: 
% Right Hand, Index
%   cross product Intermediate Phalange with Proximal Phalange to make a normal
%   pointing to the right (if the hand is palm down): y-axis of bone basis
jointNormal1 = cross(vecFingers(2:end,:), vecFingers(1:end-1,:),2);   
jointNormal1 = jointNormal1./repmat(vecMag(jointNormal1),1,3);

jointNormal1(end+1,:) = jointNormal1(end,:);

% Remove the Distal Phalange to Proximal Phalange calculation and replaces
% it with a copy of the Intermediate Phalange to Distal Phalange
% calculation. All three joints in a finger will have the same first normal
% due to the constraints in rotation of the fingers (No finger roll).
cullIdx = ~mod(1:length(jointNormal1),3);
jointNormal1(cullIdx, :) = jointNormal1(fliplr(cullIdx), :);

% Index of negative flexions (extensions)
negIdx = ~(sign(dot(jointNormal1, repmat(handy,length(jointNormal1),1), 2))+1);

jointNormal1(negIdx,:) = jointNormal1(negIdx,:) .*-1;

% Create the second normal (x-axis bone basis) for each joint.
jointNormal2 = cross(vecFingers, jointNormal1,2);
jointNormal2 = jointNormal2./repmat(vecMag(jointNormal2),1,3);

% Construct Rotation matrix for each bone
rotMat = [permute(jointNormal2', [1,3,2]),...
          permute(jointNormal1', [1,3,2]),...
          permute(vecFingers',   [1,3,2])];

rotMatDecom = permute(rotAngles(rotMat), [3,2,1]);

% handCenteredRot = rotMatDecom + repmat(handRotDecom, length(rotMatDecom),1);
handCenteredRot = rotMatDecom + repmat([0,90,90], length(rotMatDecom),1);
fingerCenteredRot = handCenteredRot(2:end,:) - handCenteredRot(1:end-1,:);
fingerCenteredRot = fingerCenteredRot(~~mod(1:length(fingerCenteredRot),3),:);

% calcAngles2(:,[3,4]) = real(reshape(fingerCenteredRot(:,2)',2,5)');
calcAngles(:,[2,1]) = handCenteredRot(fliplr(~mod(1:length(handCenteredRot),3)),2:3).*-1;
toc

%% Plot Result

if plotOn
    figure(2)
    
    for i = 1:4
        
        if i == 1
            ThumbMarkers2(1,:)  = [-1.5, -0.5, 0] + HB;
            IndexMarkers2(1,:)  = [-1, 1, 0] + HB;
            MiddleMarkers2(1,:) = [ 0, 1, 0] + HB;
            RingMarkers2(1,:)   = [ 1, 0.95, 0] + HB;
            LittleMarkers2(1,:) = [ 2, 0.80, 0] + HB;
        else
            ThumbMarkers2(i,:)  = (rotMatz(calcAngles(1,1))*rotMatx(sum(calcAngles(1,2:i))) * [0;-1;0] .* [1;-1;1] + ThumbMarkers2(i-1,:)')';
            IndexMarkers2(i,:)  = (rotMatz(calcAngles(2,1))*rotMatx(sum(calcAngles(2,2:i))) * [0;-1;0] .* [1;-1;1] + IndexMarkers2(i-1,:)')';
            MiddleMarkers2(i,:) = (rotMatz(calcAngles(3,1))*rotMatx(sum(calcAngles(3,2:i))) * [0;-1;0] .* [1;-1;1] + MiddleMarkers2(i-1,:)')';
            RingMarkers2(i,:)   = (rotMatz(calcAngles(4,1))*rotMatx(sum(calcAngles(4,2:i))) * [0;-1;0] .* [1;-1;1] + RingMarkers2(i-1,:)')';
            LittleMarkers2(i,:) = (rotMatz(calcAngles(5,1))*rotMatx(sum(calcAngles(5,2:i))) * [0;-1;0] .* [1;-1;1] + LittleMarkers2(i-1,:)')';
            
        end % END IF
        
    end % END FOR
    
    hold on
    
    scatter3(0,0,0, 'k');
    plot3([0, ThumbMarkers2(1,1)], [0, ThumbMarkers2(1,2)], [0,ThumbMarkers2(1,3)], 'k')
    plot3([0, IndexMarkers2(1,1)], [0, IndexMarkers2(1,2)], [0,IndexMarkers2(1,3)], 'k')
    plot3([0, MiddleMarkers2(1,1)], [0, MiddleMarkers2(1,2)], [0,MiddleMarkers2(1,3)], 'k')
    plot3([0, RingMarkers2(1,1)], [0, RingMarkers2(1,2)], [0,RingMarkers2(1,3)], 'k')
    plot3([0, LittleMarkers2(1,1)], [0, LittleMarkers2(1,2)], [0,LittleMarkers2(1,3)], 'k')
    
    scatter3(ThumbMarkers2(:,1), ThumbMarkers2(:,2), ThumbMarkers2(:,3), 'r');
    plot3(ThumbMarkers2(:,1), ThumbMarkers2(:,2), ThumbMarkers2(:,3), 'r')
    
    scatter3(IndexMarkers2(:,1), IndexMarkers2(:,2), IndexMarkers2(:,3), 'g');
    plot3(IndexMarkers2(:,1), IndexMarkers2(:,2), IndexMarkers2(:,3), 'g')
    
    scatter3(MiddleMarkers2(:,1), MiddleMarkers2(:,2), MiddleMarkers2(:,3), 'b');
    plot3(MiddleMarkers2(:,1), MiddleMarkers2(:,2), MiddleMarkers2(:,3), 'b')
    
    scatter3(RingMarkers2(:,1), RingMarkers2(:,2), RingMarkers2(:,3), 'm');
    plot3(RingMarkers2(:,1), RingMarkers2(:,2), RingMarkers2(:,3), 'm')
    
    scatter3(LittleMarkers2(:,1), LittleMarkers2(:,2), LittleMarkers2(:,3), 'c');
    plot3(LittleMarkers2(:,1), LittleMarkers2(:,2), LittleMarkers2(:,3), 'c')
    hold off
end