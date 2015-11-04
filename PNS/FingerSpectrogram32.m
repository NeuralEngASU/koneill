%% Loading data
% clear; close all; clc

% FileName = 'D:\Tyler\data\20120708-105356\20120708-105356-001';
FileName = 'D:\Tyler\data\20120710-155935\20120710-155935-001';

load([FileName,'.mat']);
load([FileName,'.ns3mat'],'-mat','Header')
load([FileName,'.ns3mat'],'-mat','C1')

ParametersStruct = parseNEV_FingerPress(NEV);
BSOn = round(ParametersStruct(1,1).BaselineOnTS/15);
BSOff = round(ParametersStruct(1,1).BaselineOffTS/15);

Header.ChannelLengthSamples = length(C1);
data = zeros(Header.ChannelCount,Header.ChannelLengthSamples);
Channel = Header.ChannelID;
UnitConversion = double(Header.MaxAnlgVal)./double(Header.MaxDigVal);
Units = Header.Units;
for k=1:Header.ChannelCount
    clc, disp(k)
    load([FileName,'.ns3mat'],'-mat',['C',num2str(Channel(k))])
    eval(['data(',num2str(k),',:) = double(C',num2str(Channel(k)),')*UnitConversion(',num2str(k),');']);
    clear(['C',num2str(Channel(k))])
end


data(1:32,:) = data(1:32,:) - repmat(mean(data(1:32,:)),32,1);
data_b = data(:,BSOn:BSOff);


%% Calculating spectrogram
for k=32%1:32
    clc; disp(k)
    chan = k;
    params.Fs = 2000; % sampling frequency
    params.fpass = [0 params.Fs/4]; % frequency of interest
    params.tapers = [3 5]; % tapers
    movingwin = [0.5,0.05];
       
    V = data(chan,:)';
    V_b = data_b(chan,:)';
    
    Fs = 2000;
    N = 10;
    F120 = 120;
    F180 = 180;
    Q = 20;
    h120 = fdesign.notch(N,F120,Q,Fs);
    Hd120 = design(h120);
    h180 = fdesign.notch(N,F180,Q,Fs);
    Hd180 = design(h180);
    
    v = filtfilt(Hd120.sosMatrix,Hd120.ScaleValues,V);
    v = filtfilt(Hd180.sosMatrix,Hd180.ScaleValues,v);
    
    v_b = filtfilt(Hd120.sosMatrix,Hd120.ScaleValues,V_b);
    v_b = filtfilt(Hd180.sosMatrix,Hd180.ScaleValues,v_b);
    
    [S,T,F] = mtspecgramc(v,movingwin,params);
    [S_b,T_b,F_b] = mtspecgramc(v_b,movingwin,params);
    
    x = data(45,floor(movingwin(1)/2*params.Fs):floor(movingwin(2)*params.Fs):end-floor(movingwin(1)/2*params.Fs));
    y = data(46,floor(movingwin(1)/2*params.Fs):floor(movingwin(2)*params.Fs):end-floor(movingwin(1)/2*params.Fs));
    z = data(47,floor(movingwin(1)/2*params.Fs):floor(movingwin(2)*params.Fs):end-floor(movingwin(1)/2*params.Fs));
    
    
    % Plotting
    lb = 10;
    ub = 500;
    
    idx = (F>=lb & F<=ub);
    f = F(idx);
    
    s = log(S(:,idx))';
    s_b = log(S_b(:,idx))';
    s = s - repmat(median(s_b,2),1,size(s,2));
    s = s./repmat(iqr(s_b,2),1,size(s,2));
    s = [s(:,1),s,s(:,end)];
    s = [s(1,:);s;s(end,:)];
    s = conv2(s,ones(3,3)/9,'same');
    s = s(2:end-1,2:end-1);
    
    s_m = median(s);
    s_m = exp(s_m);
    s_m = s_m*2 + lb;
    s_m(s_m<lb) = nan;
    
    s = exp(s);
    
    figure
    imagesc(T,f,s,[-1,2])
    axis xy
    hold on
    plot(T,x/100+lb-5,'r','linewidth',2)
    plot(T,y/100+lb-5,'g','linewidth',2)
    plot(T,z/100+lb-5,'b','linewidth',2)
    plot(T,s_m,'k','linewidth',2)
    hold off
    axis tight
    axis([0,25,lb-10,ub])
    
    R1 = abs(corrcoef(x,s_m));
    R2 = abs(corrcoef(y,s_m));
    R3 = abs(corrcoef(z,s_m));
    title([R1(1,2),R2(1,2),R3(1,2)])
    
    r(k,1) = R1(1,2);
    r(k,2) = R2(1,2);
    r(k,3) = R3(1,2);
end


%%
R = r - min(r(:));
R = R/max(r(:));

elecMap = [16,33,1,2,3,4,34,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,35,28,29,30,31,36,32];

ind = [R(:,1)',nan,nan,nan,nan];
ind = reshape(ind(elecMap),6,6)';
mid = [R(:,2)',nan,nan,nan,nan];
mid = reshape(mid(elecMap),6,6)';
rin = [R(:,3)',nan,nan,nan,nan];
rin = reshape(rin(elecMap),6,6)';

figure;
subplot(1,3,1)
imagesc(ind,[0,1])
subplot(1,3,2)
imagesc(mid,[0,1])
subplot(1,3,3)
imagesc(rin,[0,1])

title(sprintf('range = %0.0f to %0.0f; cc = %0.3f',lb,ub,mean(r(:))))


