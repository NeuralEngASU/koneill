function [ angleMat, fingerRotMat, tipCoord ] = Points2Angle( posMat, handData )

% Mapipulate Data
tmp = permute(posMat, [2,1,3]);
posMat2 = reshape(tmp(:), [],3);

posMat3 = posMat2-repmat(handData(1,:), size(posMat2,1),1);

handRot = handData(2:4,:)';

posMat4 = zeros(size(posMat3,1), 4);

handTrans = [[handRot;0,0,0],[handData(1,:),1]'];

%%

O1 = [0,0,0]';
O2 = handData(1,:)';

M2 = handRot;
M1 = eye(3);

P1 = posMat2;

for jj = 1:length(posMat2)
    
    P2(jj,:) = handRot' * (O1 - O2 + M1*P1(jj,:)');
    
end % END FOR


%%


% tic
% for kk = 1:size(posMat3,1)
%     
%     posMat4(kk,:) = (handTrans*[posMat3(kk,:),1]')';
%     
% end % END FOR
% toc
% 
% posMat4 = posMat4(:,1:3) - repmat(handData(1,:), length(posMat4),1);


% tic
% for kk = 1:size(posMat3,1)
%     
%     posMat4(kk,:) = (handRot*posMat3(kk,:)')';
%     
% end % END FOR
% toc
% global2local = @(x) (handRot*x(:))';
% 
% posMat4 = arrayfun(@(x,y,z) global2local([x,y,z]), posMat3(:,1), posMat3(:,2) ,posMat3(:,3) , 'UniformOutput', false)';
% posMat5 = posMat4{:}
% 
% cellfun(@(x) x', posMat4, 'UniformOutput', false)

% posMat4 = 
% posMat3 = global2localcoord((posMat2-repmat(handData(1,:), size(posMat2,1),1))', 'rr', [handData(2,:), handData(3,:), handData(4,:)]);

% Angle Matrix
% ABD, PIP, IIP, DIP
angleMat = zeros(size(posMat,1), size(posMat,2));

% Create functions
rotAngles = @(x) [(atan2(x(3,2,:), x(3,3,:)) .* (180/pi)),...
                  -acos(x(3,1,:)).*(180/pi),...  %(atan2(-x(3,1,:), sqrt(x(3,2,:).^2 + x(3,3,:).^2)) .* (180/pi)),...
                  (atan2(x(2,1,:), x(1,1,:)) .* (180/pi))];

% Find magnitude of vectors. n-dimensional norm()
vecMag = @(x) sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2);
              
%% Begin Analysis
% Finger Flexion Joints (IIP, DIP)

% Calculate the vectors from one joint to another
vecFingers = P2(2:end,:) - P2(1:end-1,:);
vecFingers(~mod(1:length(vecFingers), 4)',:) = []; % Cull excess vectors (eg: IndexTip to Middle Knuckle)

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
angleMat(:,3:4) = real(reshape(anglesTmp',2,5)');


%% Knuckle Analysis

% Idx - Mid, Idx-Rng
handy = handData(3,:);
handy = handy./norm(handy);

% Project Hand Center
% basisPlane=null(handy); %basis for the plane
% basisCoefficients= bsxfun(@minus,palmCenter,posMat2(5,:))*basisPlane ;
% handProj=bsxfun(@plus,basisCoefficients*basisPlane.', posMat2(5,:));

handz = handData(4,:);
handz = handz./norm(handz);

% handx
handx = handData(2,:);
handx = handx./norm(handx);

% handRot = [handx', handy', handz'];

handRotDecom = rotAngles(handRot);

% angleMat = handRot;

%% Plot basis
% figure(2)
% 
% tmpData = permute(P2, [1,3,2]); 
% 
% tmpData2 = reshape(tmpData(:,1,1), 4,5)';
% tmpData2(:,:,2) = reshape(tmpData(:,1,2), 4,5)';
% tmpData2(:,:,3) = reshape(tmpData(:,1,3), 4,5)';
% 
% plot3(tmpData2(1,:,1), tmpData2(1,:,2), tmpData2(1,:,3), '.k')
% hold on
% plot3(tmpData2(1,:,1), tmpData2(1,:,2), tmpData2(1,:,3), 'k')
% plot3(tmpData2(2,:,1), tmpData2(2,:,2), tmpData2(2,:,3), '.r')
% plot3(tmpData2(2,:,1), tmpData2(2,:,2), tmpData2(2,:,3), 'r')
% plot3(tmpData2(3,:,1), tmpData2(3,:,2), tmpData2(3,:,3), '.g')
% plot3(tmpData2(3,:,1), tmpData2(3,:,2), tmpData2(3,:,3), 'g')
% plot3(tmpData2(4,:,1), tmpData2(4,:,2), tmpData2(4,:,3), '.b')
% plot3(tmpData2(4,:,1), tmpData2(4,:,2), tmpData2(4,:,3), 'b')
% plot3(tmpData2(5,:,1), tmpData2(5,:,2), tmpData2(5,:,3), '.m')
% plot3(tmpData2(5,:,1), tmpData2(5,:,2), tmpData2(5,:,3), 'm')
% 
% 
% 
% 
% plot3([0, fingerx(2,1)*3] + tmpData2(2,1,1), [0, fingerx(2,2)*3]+tmpData2(2,1,2), [0, fingerx(2,3)*3]+tmpData2(2,1,3), 'r')
% hold on
% plot3([0, fingery(2,1)*3] + tmpData2(2,1,1), [0, fingery(2,2)*3]+tmpData2(2,1,2), [0, fingery(2,3)*3]+tmpData2(2,1,3), 'g')
% plot3([0, fingerz(2,1)*3] + tmpData2(2,1,1), [0, fingerz(2,2)*3]+tmpData2(2,1,2), [0, fingerz(2,3)*3]+tmpData2(2,1,3), 'b')
% 
% plot3(0,0,0,'.c')
% 
% % tmpHandData(1,:) = (handRot * handData(2,:)')';
% % tmpHandData(2,:) = (handRot * handData(3,:)')';
% % tmpHandData(3,:) = (handRot * handData(4,:)')';
% 
% tmpHandData(1,:) = handRot' * (O1-O2+M1*(handData(2,:)+handData(1,:))');
% tmpHandData(2,:) = handRot' * (O1-O2+M1*(handData(3,:)+handData(1,:))');
% tmpHandData(3,:) = handRot' * (O1-O2+M1*(handData(4,:)+handData(1,:))');
% 
% plot3([0, tmpHandData(1,1)*3], [0, tmpHandData(1,2)*3], [0, tmpHandData(1,3)*3], 'r')
% hold on
% plot3([0, tmpHandData(2,1)*3], [0, tmpHandData(2,2)*3], [0, tmpHandData(2,3)*3], 'g')
% plot3([0, tmpHandData(3,1)*3], [0, tmpHandData(3,2)*3], [0, tmpHandData(3,3)*3], 'b')
% 
% 
% dist = 100;
% xlim([tmpData2(2,1,1)-dist,tmpData2(2,1,1)+dist])
% ylim([tmpData2(2,1,2)-dist,tmpData2(2,1,2)+dist])
% zlim([tmpData2(2,1,3)-dist,tmpData2(2,1,3)+dist])
% 
% % plot3([posMat4(5,1),posMat4(6,1)],[posMat4(5,2),posMat4(6,2)],[posMat4(5,3),posMat4(6,3)],'k')
% 
% hold off
%%



% Create Normals between bones. Ex: 
% Right Hand, Index
%   cross product Intermediate Phalange with Proximal Phalange to make a normal
%   pointing to the right (if the hand is palm down): x-axis of bone basis
% % jointNormal1 = cross(vecFingers(2:end,:), vecFingers(1:end-1,:),2);   
% % jointNormal1 = jointNormal1./repmat(vecMag(jointNormal1),1,3);
% % 
% % jointNormal1(end+1,:) = jointNormal1(end,:);

%
fingerx(2,:) = cross(P2(6,1:3) - P2(5,1:3), P2(6,1:3)-P2(7,1:3));
fingerx(2,:) = fingerx(2,:)./vecMag(fingerx(2,:));

fingerx = cross(vecFingers(2:end,:), vecFingers(1:end-1,:),2);  
fingerx = fingerx./repmat(vecMag(fingerx),1,3);

fingerx(end+1,:) = fingerx(end,:);

% Remove the Distal Phalange to Proximal Phalange calculation and replaces
% it with a copy of the Intermediate Phalange to Distal Phalange
% calculation. All three joints in a finger will have the same first normal
% due to the constraints in rotation of the fingers (No finger roll).
% % cullIdx = ~mod(1:length(jointNormal1),3);
% % jointNormal1(cullIdx, :) = jointNormal1(fliplr(cullIdx), :);

cullIdx = ~mod(1:length(fingerx),3);
fingerx= fingerx(fliplr(cullIdx), :);

tmpZ = vecFingers./repmat(vecMag(vecFingers),1,3);

fingerz = -tmpZ(fliplr(~mod(1:length(vecFingers),3)),:);

fingery = cross(fingerz, fingerx);
%%

% Index of negative flexions (extensions)
% negIdx = ~(sign(dot(jointNormal1, repmat(handy,length(jointNormal1),1), 2))+1);

% Construct Rotation matrix for each bone
fingerRotMat = [permute(fingerx', [1,3,2]),...
                permute(fingery', [1,3,2]),...
                permute(fingerz', [1,3,2])];

rotMatDecom = permute(rotAngles(fingerRotMat), [3,2,1]);

% handCenteredRot = rotMatDecom + repmat(handRotDecom, length(rotMatDecom),1);
handCenteredRot = rotMatDecom + repmat([0,90,0], length(rotMatDecom),1);
handCenteredRot(:,1) = handCenteredRot(:,1) * -1;
% handCenteredRot(:,2) = (handCenteredRot(:,2));
% fingerCenteredRot = handCenteredRot(2:end,:) - handCenteredRot(1:end-1,:);
% fingerCenteredRot = fingerCenteredRot(~~mod(1:length(fingerCenteredRot),3),:);

% calcAngles2(:,[3,4]) = real(reshape(fingerCenteredRot(:,2)',2,5)');
angleMat(:,[2,1]) = handCenteredRot(:,1:2);

tipCoord = P2(4:4:end,:) - P2(1:4:end,:);

end % END FUNCTION

% EOF