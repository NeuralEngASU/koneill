% Load saved files
load('CSCData_YYYY-MM-DD_hh-mm-ss_CSC.mat')
load('SpikeData_YYYY-MM-DD_hh-mm-ss_SE.mat')

CSCData1 = CSCData;
SpikeData1 = SpikeData;

load('CSCData_YYYY-MM-DD_hh-mm-ss_CSC.mat')
load('SpikeData_YYYY-MM-DD_hh-mm-ss_SE.mat')

CSCData2 = CSCData;
SpikeData2 = SpikeData;

% Find trial events
[ trialEvents1, ~ ] = FindLaserEvents( CSCData1 );

% RasterComp
RasterComp(SpikeData1, SpikeData2, saveData, [5,7,9], 1, trialEvents1, 'b', trialEvents1(:,1), [-5, 100] , CSCData1.freq )

