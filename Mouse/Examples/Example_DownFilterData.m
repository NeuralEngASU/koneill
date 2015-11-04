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

TrialRaw( CSCData, 1, 8, 0, [0,1], 'k', [5*60, 10*60], [])

