%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CalcPLI2
%   Author: Kevin O'Neill
%   Date: 2015.03.16
%
%   Desc: 
%           Compute the Phase-Lag Index (PLI) and confidence between each 
%           pair of channels. The PLI is defined as the absolute value of
%           the mean sign of the instantaneous phase difference between a
%           pair of signals:
%
%                   PLI = |<sign(deltaPhi(t))>|
%
%                   PLI = abs(mean(deltaPhi));
%
%                           range = [0,1]
%
%           PLI measures the asymmetry in the distribution of instantaneous
%           phase differences between two signals. If the signals are
%           coupled then the PLI will return a value close to 1 ( or PLI
%           will be 0 if the signals are coupled around 0 mod pi. In these
%           cases the signals are 0, pi, 2pi, ... radians out of phase).
%
%           If signals are not coupled, then PLI will return a value close
%           to 0. If signals are not coupled this implies that there is
%           either no, or a very weak common signal (volume conduction)
%           that influences both signals.
%
%           Surrogate data (randomly phase shifted original data) are used
%           to test for confidence. A z-score calculation is sufficient to
%           check if the calculated PLI is signifigant. The sign of the
%           z-score tells us if the signals are significantly coupled or
%           non-coupled.
%           
%           The output of PLI for two sine waves with phase differences
%           from 0:pi is approximated with the code below:
%
<<<<<<< HEAD
%           tIdx = 1:1000;
%           phaseAngle = linspace(0, pi, length(tIdx));
%           sig1 = sin(pi*tIdx/100);
%           PLI = zeros(1,length(tIdx));
%           for kk = 1:length(phaseAngle)
%                 sig2 = sin(pi*tIdx/100 + phaseAngle(kk));
% 
%                 sig1Phi = atan(imag(hilbert(sig1))./sig1);
%                 sig2Phi = atan(imag(hilbert(sig2))./sig2);
% 
%                 deltaPhi = sig1Phi - sig2Phi;
%                 PLI(kk) = abs(mean(sign(deltaPhi)));
%           end
%           plot(phaseAngle, PLI); title('PLI v phase'); xlabel('Phase, radians'); ylabel('PLI');
=======
figure(2)
          tIdx = 1:1000;
          phaseAngle = linspace(0, pi, length(tIdx));
          sig1 = sin(pi*tIdx/100);
          PLI = zeros(1,length(tIdx));
          for kk = 1:length(phaseAngle)
                sig2 = sin(pi*tIdx/100 + phaseAngle(kk));

                sig1Phi = atan(imag(hilbert(sig1))./sig1);
                sig2Phi = atan(imag(hilbert(sig2))./sig2);

                deltaPhi = sig1Phi - sig2Phi;
                PLI(kk) = abs(mean(sign(deltaPhi)));
          end
          plot(phaseAngle, PLI); title('PLI v phase'); xlabel('Phase, radians'); ylabel('PLI');
>>>>>>> origin/master
%
%
%       Usage:
%           [R, PLI, zScore, phaseShift] = CalcPLI2(pathName,NEVFileName,tStart,tSpan);
%
%       Outputs:
%           Returns the Phase Similarity, PLI, Z-Score, and Phase Shift
%           matricies to the workspace and saves a .mat file to the NEV
%           file location.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R, PLI, zScore, phaseShift] = CalcPLI2(pathName,NEVFileName, tStart, tSpan)

%% Load Data %%

% Make sure the data files are in MATLAB's path.

expr = '(nev)$'; % Expresstion for regex

% Make sure the file given is for a .nev file
testName = regexp(NEVFileName, expr, 'Tokens');

fullNEVPath = fullfile(pathName, NEVFileName);

% Test if the file is a .nev
if isempty(testName)
    warning('File name is incorrect. Please send the full file path to the NEV file.')
end % END IF testName

% Test to make sure the file exists
if ~(exist(fullNEVPath, 'file') == 2)
    warning('File does not exist. Please send the full file path to the NEV file.')
end % END IF testName

% Load the data from the NEV and nsx files
disp('##### Loading Data #####')

NEV = openNEV(fullNEVPath);

if exist([fullNEVPath(1:end-3),'ns4'],'file')
    
    nsx2mat([fullNEVPath(1:end-3),'ns4'])
    load([fullNEVPath(1:end-3),'ns4mat'], '-mat')
    
elseif exist([fullNEVPath(1:end-3),'ns5'],'file')
        
    nsx2mat([fullNEVPath(1:end-3),'ns5'])
    load([fullNEVPath(1:end-3),'ns5mat'], '-mat')
        
end % END IF

disp('##### Data Loaded #####')
%% Make files names for saved data

fullNameMat = fullfile(pathName, [NEV.MetaTags.Filename, '_PLI.mat']);

%% NS5 Output

outputStruct = openNSx([fullNEVPath(1:end-3),'ns5'], 'read', 'c:1:32','t:15000000:44999999');


for aa = 1:32
    
    eval(['C', num2str(aa), ' = outputStruct.Data(', num2str(aa), ',:);']);
    
end % END ns5

%% Phase Lag Index

if tStart == 0
    idxStart = 1;
    idxSpan  = tSpan * Header.Fs;
    
    idxPLI   = idxStart:(idxSpan+idxStart-1);
else
    idxStart = tStart * Header.Fs;
    idxSpan  = tSpan  * Header.Fs;
    
    idxPLI   = idxStart:(idxSpan+idxStart-1);
end

% Allocate variables
sig1 = zeros(1,length(idxPLI)); % Signal 1 container
sig2 = zeros(1,length(idxPLI)); % Signal 2 container

numChannels  = 32; %Header.ChannelCount; % Number of channels
spanChannels = 1:numChannels; % The total channel span
spanComp     = spanChannels; % The span of the channel comparisons

surrDataCount = 250; % Number of surrogate data trials

surrData = zeros(1,length(sig1)); % Surrogate data container
surrFFT  = zeros(1,length(sig1)); % Surrogate FFT container
sig1Phi  = zeros(1,length(sig1)); % Signal 1 phase container
sig2Phi  = zeros(1,length(sig1)); % Signal 2 phase container

deltaPhi = zeros(1,length(sig1)); % Phase difference container

% Randomly generated phase shifts. Only use the upper half of the matrix
phaseShift = -pi + (pi+pi)*rand(numChannels,numChannels,surrDataCount);

% Result matricies
R      = zeros(numChannels, numChannels, surrDataCount+1); % Phase similarity container
PLI    = zeros(numChannels, numChannels, surrDataCount+1); % Phase-lag index container
zScore = zeros(numChannels, numChannels, surrDataCount+1); % Z-Score container

% Randians to samples conversion
T = 1/Header.Fs;
radSamp = floor(phaseShift./(2*pi) * ((2*pi)/T));

timeTaken = 0;

parpool(10);

tSet = tic; % Set a overall start time
for ii = spanChannels % Loop over each channel
    
    % Copy the raw data into new variables
%     eval(['sig1 = double(C', num2str(ii),'(idxPLI));']);
    eval(['sig1 = double(C', num2str(ii),');']); % NS5
    
    % Calculate the instantaneous phase of signal 1
    sig1Phi = atan(imag(hilbert(sig1))./sig1);
        
    for jj = spanComp % Loop over all comparisons (1:32, 2:32, 3:32, ..., 31:32, 32)
        
        % Copy the raw data into new variables
%         eval(['sig2 = double(C', num2str(jj),'(idxPLI));']);
        eval(['sig2 = double(C', num2str(jj),');']); % NS5
        
        % Calculate the instantaneous phase of signal 2
        sig2Phi = atan(imag(hilbert(sig2))./sig2);
        
        % Start loop time
        tLoop = tic; % Start a loop time counter
        
        parfor kk = 1:surrDataCount % loop over the surrogate data trial

            % Shift the phase of the surrogate data 
            surrPhi = circshift(sig2Phi', radSamp(ii,jj,kk));
            
            % Calculate the instantaneous phase difference between signal 1
            % and the surrogate signal 2
            deltaPhi = sig1Phi - surrPhi';
            
            % Save the phase-similarity and phase-lag index scores to their
            % matricies
            R(ii,jj,kk) = abs(1/length(deltaPhi) * sum( exp(1i * deltaPhi(1:end-1))) );
            PLI(ii,jj,kk) = abs(mean(sign(deltaPhi)));
            
        end % END FOR surrDataCount
               
        % Calculate the instantaneous phase difference between signal 1 and
        % signal 2
        deltaPhi = sig1Phi - sig2Phi;
        
        % Compute and store the phase similarity and the phase-lag index
        % in their own matrix. This computation is for the 'Real' signals
        R(ii,jj,end) = abs(1/length(deltaPhi) * sum( exp(1i .* deltaPhi)) );
        PLI(ii,jj,end) = abs(mean(sign(deltaPhi)));
        
        % Compute the zScore for each PLI. A |score| >= 1.96 is above the
        % 95% confidence interval. 
        zScore(ii,jj,:) = zscore(PLI(ii,jj,:) - PLI(ii,jj,end));
        
        %%%% User Outputs %%%%
        timeSpent = toc(tLoop); % Capture the time spent computing the PLI
        
        % Calculate the total number of computations that have
        % completed
        if ii == 1
            numComplete = sum(spanComp<jj)*surrDataCount + surrDataCount;
        elseif ii > 1
            numComplete = sum(spanChannels(end - (ii-2):end))*surrDataCount + sum(spanComp<jj)*surrDataCount + surrDataCount;
        end % END IF ii
        
        numTotal =  sum(1:32)*surrDataCount; % Calculate the toal number of loops that will occur
        
        timeLeft = (numTotal - numComplete)/surrDataCount * timeSpent; % Approximation of the time left to finish code
        timeTaken = timeTaken + timeSpent; % Approximation of the time spent running the code so far
        
        % Display user readable outputs about how many minutes are left
        clc;
        disp(sprintf('%d/%d', numComplete, numTotal));
        disp(sprintf('Time left: %f, minutes',timeLeft/60))
        disp(sprintf('Time spent: %f, minutes',timeTaken/60))
        disp(fullNEVPath)
        %%%% END User Outputs %%%%
        
        
    end % END FOR sig2 spanComp
    
    % Remove the first comparison after each cycle. This allows us to
    % compute only the upper half of the confusion matrix.
    spanComp(1) = [];
    
end % END FOR sig1 spanChannels

% Flip the third dimension in order to bring the 'real' value to the first
% plane
R      = R(:,:,end:-1:1);
PLI    = PLI(:,:,end:-1:1);
zScore = zScore(:,:,end:-1:1);

% Save Data
save(fullNameMat, 'R', 'PLI', 'zScore', 'phaseShift');


poolobj = gcp('nocreate');
delete(poolobj);

% Display the total time required to run the code.
disp(sprintf('Total time to complete: %d', toc(tSet)/60));
end % END FUNCTION
% EOF