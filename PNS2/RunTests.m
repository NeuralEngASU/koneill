%% Load in Data File

%% Down Filter

% Downsample and filter for Spike data. Spikes are around 1000 Hz
DownFilter( 1, TDTData, 1, 2000, 1, 'highpass', 500);

%% Plot Trial Raw

PNSTrialRaw(); % Plot a small section of a channel to make sure the data looks good.

%% Find Events

events = PNSFindEvents(); % Looks for event headers. Index flex, middle flex, ect...

events = events ./ 1; % Make sure that the timestamps between the downsampled data and events match.

eventNames = {{'ThumbFlex'},...
              {'IndexFlex'},...
              {'MiddleFlex'},...
              {'RingFlex'},...
              {'LittleFlex'},...
              {'ThumbExt'},...
              {'IndexExt'},...
              {'MiddleExt'},...
              {'RingExt'},...
              {'LittleExt'},...
              {'ThumbAbd'},...
              {'IndexAbd'},...
              {'MiddleAbd'},...
              {'RingAbd'},...
              {'LittleAbd'}};

%% Plot Raw Data for Each movement

% Plot finger movement paths:
fingerMove = true;

for i = 1:13
    PNSTrialRaw( fingerMove,['Name_',eventNames{i}]);
end % END FOR

%% Plot Raster data for each Finger movement

for i = 1:13
    PNSRastergram( events(:,i),['Name_',eventNames{i}]);
end % END FOR

%% Plot PSTH data

for i = 1:13
    PNSPSTH( events(:,i),['Name_',eventNames{i}]);
end % END FOR

%% Plot cross correllation for each finger movement and print save the top correlated electrodes

for i = 1:13
    R(:,:,i) = PNSCorrCoef( 1, events(:,i),['Name_',eventNames{i}]);
    PNSCorrCoef( 2, events(:,i),['Name_',eventNames{i}]); % Plot both plot options
end % END FOR

save('CrossCorrellations.mat',R);

%% Plot Spike Data

electrodes = 1:96;
channels = 1:96;
for i = channels
    PNSSpikeWaveform(i,1);
    PNSSpikeWaveform(i,2);
    SNR(i) = PNSSNR(i,baseline);
end % END FOR

%% Plot spectrum

%% Plot heatmap

for i = 1:13
    GenHeatmap(data(events(:,i)), ['Name_',eventNames{i}]);
    
end % END FOR