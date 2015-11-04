function [xhat,P] = kalman_test(xhat,P,z,TRAIN)

% subtract mean during training
% z = z-TRAIN.zmu;

% step 1: time-update equations
xhatm = TRAIN.A*xhat(:,1); %previous xhat
Pm = TRAIN.A*P(:,:,1)*TRAIN.A' + TRAIN.W; %previous P

% step 2: measurement-update equations
K = Pm*TRAIN.H'*pinv(TRAIN.H*Pm*TRAIN.H'+TRAIN.Q);
xhat(:,2) = xhatm + K*(z-TRAIN.H*xhatm); %current xhat
P(:,:,2) = (eye(size(K,1))-K*TRAIN.H)*Pm; %current P

% rotate buffers (discard old data)
P(:,:,1) = P(:,:,2);
xhat(:,1) = xhat(:,2);