%% DownFilter() and TrialRaw()
% Pull data from the saved NeuraLynx files
nsc2mat(1);

% Load saved file
load('CSCData_YYYY-MM-DD_hh-mm-ss_CSC.mat')

% Down sample for LFP data. LFPs are around <400 Hz
DownFilter( 1, CSCData, 1, 1000, 0, '', 0)


% Down sample and filter for Spike data. Spikes are around 1000 Hz
DownFilter( 1, CSCData, 1, 2000, 1, 'highpass', 500)

% Plot spike raw data
load('CSCData_YYYY-MM-DD_hh-mm-ss_CSC_ds2000_fsH500.mat')

TrialRaw( CSCData, 1, 8, 0, [50, 60], 'k', [10, 110], [])

%% MousePSTH()

% Load data file
load('SpikeData_YYYY-MM-DD_hh-mm-ss_SE.mat')

% Run PSTH for a specified time window (5 to 10 seconds)
MousePSTH( SpikeData, 1, 0, [], '', [5, 10], [], 32000)

% Run PSTH for a specified window around several events (-20 to +20 around
% event)
MousePSTH( SpikeData, 1, 0, [], '', [-20, 20], [100, 200,300, 400], 32000)

%% SpikeWaveform()

% Plot the spike waveform for a selected channel
SpikeWaveform( SpikeData, 1, 1, 0); % Channel 1
SpikeWaveform( SpikeData, 1, 2, 0); % Channel 2

%% MouseSNR()

% Find the Signal-to-Noise Ratio for a selected channel. This function uses
% a baseline to measure noise (from 1-10 seconds)
MouseSNR( CSCData, SpikeData, 1, [1, 10] ); % Channel 1

MouseSNR( CSCData, SpikeData, 2, [1, 10] ); % Channel 2

%% Rastergram

% Plots a rastergram from channel 1, from 10 seconds to 110 seconds
Rastergram( SpikeData, 1, 1, 0, [], '', [10, 110], [], 0 )

% Plots a rastergram from channel 2, for -20 to +20 seconds around events
Rastergram( SpikeData, 1, 2, 0, [], '', [-20, +20], [100, 200, 300, 400], 0 )

% Plots a rastergram from channel 2, for -20 to +20 seconds around events
% and displays a blue trial identifier (from -10 to 0) in the background
Rastergram( SpikeData, 1, 2, 1, [-10, 0], 'b', [-20, +20], [100, 200, 300, 400], 0 )


%% MouseCorrCoef()

% Outputs the complete cross correlation coefficient matric for the
% specified time segment.
R = MouseCorrCoef( CSCData, 1, 0, [10, 110], [], '100SecCorrCoef', 1, 2);

% Outputs the complete cross correlation coefficient matric for the
% specified time segment, including the laser data.
R = MouseCorrCoef( CSCData, 1, 1, [10, 110], [], '100SecCorrCoefLaser', 1, 2);

% Outputs the complete cross correlation coefficient matric for the
% specified time segment, including the laser data, around events. In the
% figure the R value will be mean'd
R = MouseCorrCoef( CSCData, 1, 1, [-20, 20], [100, 200, 300, 400], '100SecCorrCoefLaserEvents', 2, 2);
