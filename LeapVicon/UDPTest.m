%%

hudpr = dsp.UDPReceiver('LocalIPPort', 9001, 'MaximumMessageLength', 1024);

data = {};
% for kk = 1:100
kk = 1;
while 1
    
    data{kk} =  step(hudpr);
    kk = kk+1;
    pause(1/60);
end % END FOR

release(hudpr);

%% Regexp

jointDataRight = zeros(5,4,3, length(data)); % Position (finger by joint (proximal, intermediate, distal) by x-y-z
jointDataLeft = zeros(5,4,3,length(data));

handRight = zeros(4,3, length(data)); % PalmPosisiton; xbasis; ybasis; zbasis 
handLeft = zeros(1,3, length(data)); % Position 

wristRight = zeros(1,3);
wristLeft = zeros(1,3);


angleMatOut = zeros(5, 4, length(data));

for ii = 2:length(data)/2
    
    tmpString = char(data{ii}');
    
    if ~strcmp(tmpString, 'NoHand')
    
        C = tmpString(2:end);
        eval(C)
        
        if exist('RightHandJointPos')
            tic
            handRight(:,:,ii) = RightHandOrn(:,[1,3,2]);
            handRight(:,2,ii) = -1*handRight(:,2,ii);
            
            
            jointDataRight(:,:,1, ii) = reshape(RightHandJointPos(1:3:end), 4, 5)';
            jointDataRight(:,:,2, ii) = -1*reshape(RightHandJointPos(3:3:end), 4, 5)';
            jointDataRight(:,:,3, ii) = reshape(RightHandJointPos(2:3:end), 4, 5)';
            
            [angleMatOut(:,:,ii), fingerRot(:,:,:,ii), tipCoord(:,:,ii)] = Points2Angle(jointDataRight(:,:,:,ii), handRight(:,:,ii));
                        
%             plot3(jointDataRight(1,:,1,ii), jointDataRight(1,:,2,ii), jointDataRight(1,:,3,ii), '.k')
%             hold on
%             plot3(jointDataRight(1,:,1,ii), jointDataRight(1,:,2,ii), jointDataRight(1,:,3,ii), 'k')
%             plot3(jointDataRight(2,:,1,ii), jointDataRight(2,:,2,ii), jointDataRight(2,:,3,ii), '.r')
%             plot3(jointDataRight(2,:,1,ii), jointDataRight(2,:,2,ii), jointDataRight(2,:,3,ii), 'r')
%             plot3(jointDataRight(3,:,1,ii), jointDataRight(3,:,2,ii), jointDataRight(3,:,3,ii), '.g')
%             plot3(jointDataRight(3,:,1,ii), jointDataRight(3,:,2,ii), jointDataRight(3,:,3,ii), 'g')
%             plot3(jointDataRight(4,:,1,ii), jointDataRight(4,:,2,ii), jointDataRight(4,:,3,ii), '.b')
%             plot3(jointDataRight(4,:,1,ii), jointDataRight(4,:,2,ii), jointDataRight(4,:,3,ii), 'b')
%             plot3(jointDataRight(5,:,1,ii), jointDataRight(5,:,2,ii), jointDataRight(5,:,3,ii), '.m')
%             plot3(jointDataRight(5,:,1,ii), jointDataRight(5,:,2,ii), jointDataRight(5,:,3,ii), 'm')
%             plot3( handRight(:,1,ii),  handRight(:,2,ii),  handRight(:,3,ii), '.c')
%             plot3(0,0,0, 'ko')
%             plot3(0,0,0, 'k.')
% %             plot3(handProj(:,1,ii), handProj(:,2,ii), handProj(:,3,ii), 'oc')
%             
% %             plot3([handProj(:,1,ii), handProj(:,1,ii)+ angleMatOut(1, 1, ii)*3], [handProj(:,2,ii), handProj(:,2,ii)+ angleMatOut(2, 1, ii)*3], [handProj(:,3,ii), handProj(:,3,ii)+ angleMatOut(3, 1, ii)*3], 'r')
% %             plot3([handProj(:,1,ii), handProj(:,1,ii)+ angleMatOut(1, 2, ii)*3], [handProj(:,2,ii), handProj(:,2,ii)+ angleMatOut(2, 2, ii)*3], [handProj(:,3,ii), handProj(:,3,ii)+ angleMatOut(3, 2, ii)*3], 'g')
% %             plot3([handProj(:,1,ii), handProj(:,1,ii)+ angleMatOut(1, 3, ii)*3], [handProj(:,2,ii), handProj(:,2,ii)+ angleMatOut(2, 3, ii)*3], [handProj(:,3,ii), handProj(:,3,ii)+ angleMatOut(3, 3, ii)*3], 'b')
% 
%             tmpHandRight(:,:,ii) = handRight(:,:,ii)*3 + repmat(handRight(1,:,ii),4,1);
%             plot3([handRight(1,1,ii), tmpHandRight(2,1,ii)], [handRight(1,2,ii), tmpHandRight(2,2,ii)], [handRight(1,3,ii), tmpHandRight(2,3,ii)], 'r')
%             plot3([handRight(1,1,ii), tmpHandRight(3,1,ii)], [handRight(1,2,ii), tmpHandRight(3,2,ii)], [handRight(1,3,ii), tmpHandRight(3,3,ii)], 'g')   
%             plot3([handRight(1,1,ii), tmpHandRight(4,1,ii)], [handRight(1,2,ii), tmpHandRight(4,2,ii)], [handRight(1,3,ii), tmpHandRight(4,3,ii)], 'b')   
% 
%             xlim([-250,250])
%             ylim([-250,250])
%             zlim([0,250])
%             hold off
            
            %%%% Circle Plot %%%%
            
            O1 = [0,0,0]';
            O2 = [0,0,0]';
            
            M2 = basisOffsets;
            M1 = eye(3);
            
            P1 = tipCoord(:,:,ii);
            
            for jj = 1:5
                
                tipDist(jj,:) = M2(:,:,jj)' * (O1 - O2 + M1*P1(jj,:)');
                
            end % END FOR


            plot(tipDist(2,1), tipDist(2,2), 'ko')
            hold on
            
            circx = sind(angleMatOut(2,1,ii));
            circy = sind(-angleMatOut(2,2,ii));
            plot(circx*75, circy*75, 'or')
            hold off
     
            
            xlim([-50, 50])
            ylim([-50, 50])
            drawnow
            timeOut = toc;
            pause(1/60-timeOut)
            
        end % END IF right hand exists   
    end % END IF hands exist    
end % END FOR time

%% Circle Plot - Position Setup

zeroIdx = sum(squeeze(tipCoord(2,1,:)) == 0) + 1;

tipOffsets = mean(tipCoord(:,:,57:end),3);
basisOffsets = mean(fingerRot(:,:,:,57:end),4);

O1 = [0,0,0]';
O2 = [0,0,0]';

M2 = basisOffsets;
M1 = eye(3);

P1 = tipOffsets;

for jj = 1:5
    
    tipDist(jj,:) = M2(:,:,jj)' * (O1 - O2 + M1*P1(jj,:)');
    
end % END FOR
% 
% %% Circle Plot - Position
% 
% 
% 
% 
% 
% 
% 
% 
% %% Circle Plot - Rotations
% 
% circx = sind(angleMatOut(2,1,:) - averageAbd);
% circy = sind(angleMatOut(2,2,:) - averageFlex);
% figure(3)
% for aa = 1:length(circx)
%     
%     plot(circx(aa), circy(aa), 'ok')
%         xlim([-1,1])
%     ylim([-1,1])
%     drawnow
%     pause(1/60)
% 
%     
% end % END FOR
% 

%%
hudpr = dsp.UDPReceiver('LocalIPPort', 9001, 'MaximumMessageLength', 1024);

cosData = cos(0:0.1:2*pi) * 10;
sinData = sin(0:0.1:2*pi)*10;
data = {};
% for kk = 1:100
ii = 1;
while 1
    tic
    data{ii} =  step(hudpr);
    
    
    
    tmpString = char(data{ii}');
    
    if ~strcmp(tmpString, 'NoHand')
    
        C = tmpString(2:end);
        eval(C)
        
        if exist('RightHandJointPos')
            
            handRight(:,:,ii) = RightHandOrn(:,[1,3,2]);
            handRight(:,2,ii) = -1*handRight(:,2,ii);
            
            
            jointDataRight(:,:,1, ii) = reshape(RightHandJointPos(1:3:end), 4, 5)';
            jointDataRight(:,:,2, ii) = -1*reshape(RightHandJointPos(3:3:end), 4, 5)';
            jointDataRight(:,:,3, ii) = reshape(RightHandJointPos(2:3:end), 4, 5)';
            
%             [angleMatOut(:,:,ii), fingerRot(:,:,:,ii), tipCoord(:,:,ii)] = Points2Angle(jointDataRight(:,:,:,ii), handRight(:,:,ii));
                
%             figure(1)
%             plot3(jointDataRight(1,:,1,ii), jointDataRight(1,:,2,ii), jointDataRight(1,:,3,ii), '.k')
%             hold on
%             plot3(jointDataRight(1,:,1,ii), jointDataRight(1,:,2,ii), jointDataRight(1,:,3,ii), 'k')
%             plot3(jointDataRight(2,:,1,ii), jointDataRight(2,:,2,ii), jointDataRight(2,:,3,ii), '.r')
%             plot3(jointDataRight(2,:,1,ii), jointDataRight(2,:,2,ii), jointDataRight(2,:,3,ii), 'r')
%             plot3(jointDataRight(3,:,1,ii), jointDataRight(3,:,2,ii), jointDataRight(3,:,3,ii), '.g')
%             plot3(jointDataRight(3,:,1,ii), jointDataRight(3,:,2,ii), jointDataRight(3,:,3,ii), 'g')
%             plot3(jointDataRight(4,:,1,ii), jointDataRight(4,:,2,ii), jointDataRight(4,:,3,ii), '.b')
%             plot3(jointDataRight(4,:,1,ii), jointDataRight(4,:,2,ii), jointDataRight(4,:,3,ii), 'b')
%             plot3(jointDataRight(5,:,1,ii), jointDataRight(5,:,2,ii), jointDataRight(5,:,3,ii), '.m')
%             plot3(jointDataRight(5,:,1,ii), jointDataRight(5,:,2,ii), jointDataRight(5,:,3,ii), 'm')
%             plot3( handRight(:,1,ii),  handRight(:,2,ii),  handRight(:,3,ii), '.c')
%             plot3(0,0,0, 'ko')
%             plot3(0,0,0, 'k.')
% %             plot3(handProj(:,1,ii), handProj(:,2,ii), handProj(:,3,ii), 'oc')
%             
% %             plot3([handProj(:,1,ii), handProj(:,1,ii)+ angleMatOut(1, 1, ii)*3], [handProj(:,2,ii), handProj(:,2,ii)+ angleMatOut(2, 1, ii)*3], [handProj(:,3,ii), handProj(:,3,ii)+ angleMatOut(3, 1, ii)*3], 'r')
% %             plot3([handProj(:,1,ii), handProj(:,1,ii)+ angleMatOut(1, 2, ii)*3], [handProj(:,2,ii), handProj(:,2,ii)+ angleMatOut(2, 2, ii)*3], [handProj(:,3,ii), handProj(:,3,ii)+ angleMatOut(3, 2, ii)*3], 'g')
% %             plot3([handProj(:,1,ii), handProj(:,1,ii)+ angleMatOut(1, 3, ii)*3], [handProj(:,2,ii), handProj(:,2,ii)+ angleMatOut(2, 3, ii)*3], [handProj(:,3,ii), handProj(:,3,ii)+ angleMatOut(3, 3, ii)*3], 'b')
% 
%             tmpHandRight(:,:,ii) = handRight(:,:,ii)*3 + repmat(handRight(1,:,ii),4,1);
%             plot3([handRight(1,1,ii), tmpHandRight(2,1,ii)], [handRight(1,2,ii), tmpHandRight(2,2,ii)], [handRight(1,3,ii), tmpHandRight(2,3,ii)], 'r')
%             plot3([handRight(1,1,ii), tmpHandRight(3,1,ii)], [handRight(1,2,ii), tmpHandRight(3,2,ii)], [handRight(1,3,ii), tmpHandRight(3,3,ii)], 'g')   
%             plot3([handRight(1,1,ii), tmpHandRight(4,1,ii)], [handRight(1,2,ii), tmpHandRight(4,2,ii)], [handRight(1,3,ii), tmpHandRight(4,3,ii)], 'b')   
% 
%             xlim([-250,250])
%             ylim([-250,250])
%             zlim([0,250])


            %%% Circle Plot %%%%
            
            O1 = [0,0,0]';
            O2 = [0,0,0]';
            
            M2 = basisOffsets;
            M1 = eye(3);
            
            P1 = tipCoord(:,:,ii);
            
            for jj = 1:5
                
                tipDist(jj,:) = M2(:,:,jj)' * (O1 - O2 + M1*P1(jj,:)');
                
            end % END FOR

            
            plot(cosData, sinData, 'Color', [.5,0.5,0.5], 'Linewidth', 2.75)
            hold on
            plot(tipDist(2,1), tipDist(2,2), '.k', 'MarkerSize', 16)
            xlim([-50,50])
            ylim([-50,50])
            axis('square')


            hold off
            drawnow
        end
    end
    
    
    timeOut = toc;
    if timeOut < (1/60)
        pause(1/60-timeOut);
    end
    ii = ii+1;
end % END FOR

release(hudpr);





