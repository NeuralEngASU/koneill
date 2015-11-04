function [PARAM,x,z]=kalman_train(x,z)

% calculate mu's, the means of kinematics and features, and subtract off
% PARAM.xmu=mean(x,2); x=x-repmat(PARAM.xmu,[1 size(x,2)]);
% PARAM.zmu=mean(z,2); z=z-repmat(PARAM.zmu,[1 size(z,2)]);

% length of the signals
M=size(x,2);

% calculate A, the state-to-state transformation (hand kinematics)
A1=x(:,2:M)*x(:,1:(M-1))';
A2=x(:,1:(M-1))*x(:,1:(M-1))';
PARAM.A=A1*pinv(A2);

% calculate W, the covariance of the noise in the kinematics
W1=x(:,2:M)*x(:,2:M)';
W2=x(:,1:(M-1))*x(:,2:M)';
PARAM.W=(1/(M-1))*(W1-PARAM.A*W2);

% cross-correlation and autocorrelations of x and z
PARAM.Pzx=z(:,1:M)*x(:,1:M)';
PARAM.Rxx=x(:,1:M)*x(:,1:M)';
PARAM.Rzz=z(:,1:M)*z(:,1:M)';

% calculate H, the transformation matrix from measured features to state
PARAM.H=PARAM.Pzx*pinv(PARAM.Rxx);

% calculate Q, the covariance of noise in the measured features
PARAM.Q=(1/M)*(PARAM.Rzz-PARAM.H*PARAM.Pzx');