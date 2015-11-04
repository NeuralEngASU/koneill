function [] = RasterPlotCorrComparison( trialStruct, trialTypes, trials, channels, trialLen, spikeTimes, H, varargin )

% trialStruct
% Trial Type
% Trials
% Spike Times
% Varargin

numIn = nargin;

%%%%% INIT %%%%%
fs         = 30000; % default sampling rate
flagLabel  = 1;     % boolean to show labels
spikeWidth = 1;     % spike thickness
spikeColor = {'r', 'g', 'b'};   % spike color
trialGap   = 1.5;   % distance between trials
buffer     = 7500;  % beginning and end buffer in samples

switch numIn
    case 7 % Use defaults on Spike Times
        if isempty(H)
            figure;
            H=gca;
        end
    case 8 % Alters sampling rate
        figure;
        fs = varargin{1};
    case 9 % Alters labels
        figure;
        fs        = varargin{1};
        flagLabel = varargin{2};
    case 10 % Spike Width
        figure;
        fs         = varargin{1};
        flagLabel  = varargin{2};
        spikeWidth = varargin{3};
    case 11 % Spike Color
        figure;
        fs         = varargin{1};
        flagLabel  = varargin{2};
        spikeWidth = varargin{3};
        spikeColor = varargin{4};
    otherwise
        error ('Invalid Arguments');
end



%%%%% MAIN %%%%%
hold on

% Patch
patchX = [0,3,3,0];
patchY = [30,30,0,0];
pData = patch(patchX, patchY, 1);

set(pData, 'FaceColor', 'k')
set(pData, 'EdgeColor', 'none')
set(pData, 'FaceAlpha', 0.15)


for j = 1:length(trialTypes)
    
    % Parses information
    trialType    = trialTypes{j};
    rasterTrials = eval(['trialStruct.', trialType,'(trials,:)']);
    numTrials    = 1;%length(trials);
    numChans     = 1;%size(channels, 1);
    
    % Loops over the number of trials
    for i = 6%numChans
        for k = 1:3%numTrials
            
%             if j==3 && k==20 && i==17
%                 disp([k, i, j])
%             end

            % Size of each trial to plot (Samples)
            trialSize = trialLen + 2* buffer;
            startTime = rasterTrials(k,4) - buffer;
            
            trialSpikesIdx = spikeTimes{channels(i)} > startTime & spikeTimes{channels(i)} < [startTime + trialSize];
            trialSpikes = spikeTimes{channels(i)}(trialSpikesIdx);
            
            % Checks if there were spikes within the trial region
            if ~isempty(trialSpikes)
                spikeNum = length(trialSpikes);
                
                xx = ones(1, 3 * spikeNum) * nan;
                yy = ones(1, 3 * spikeNum) * nan;
                
                % Convert to seconds and organize xx in to sets of 3 points
                xx(2:3:3*spikeNum) = double((trialSpikes - startTime)).*(1/fs);
                xx(1:3:3*spikeNum) = xx(2:3:3*spikeNum) - eps(1);
                xx(3:3:3*spikeNum) = xx(2:3:3*spikeNum) + eps(1);
                
                yy(2:3:3*spikeNum) = (40 - j*10) - trialGap;
                yy(1:3:3*spikeNum) = (30 - j*10);
                
                plot(H, (xx - (buffer/fs)), yy, spikeColor{j}, 'linewidth', spikeWidth)
                
            end % END IF
        end % END FOR
    end % END FOR
    
    
%     % Plots the trial start and end lines
%     xFlag = ones(1, 6*numTrials) * nan;
%     yFlag = ones(1, 6*numTrials) * nan;
%     
%     xFlag(1:6:6*numTrials) = ([1:numTrials]-1) .* (trialSize) ./fs;
%     xFlag(4:6:6*numTrials) = ([1:numTrials]) .* (trialSize) ./fs - 2*buffer/fs;
%     xFlag(2:6:6*numTrials) = xFlag(1:6:6*numTrials) + eps(1);
%     xFlag(5:6:6*numTrials) = xFlag(4:6:6*numTrials) + eps(1);
%     
%     yFlag(1:6:6*numTrials) = 0;
%     yFlag(4:6:6*numTrials) = 0;
%     yFlag(2:6:6*numTrials) = numChans*10;
%     yFlag(5:6:6*numTrials) = numChans*10;
%     
%     plot(H, xFlag + (j-1)*trialSize/fs*20, yFlag, [spikeColor{j}, ':'])
      
end % END FOR

hold off


% Styles plot
xlim([-0.5, 3.5])
ylim([0,30])

if flagLabel
    ylabel(sprintf('Finger\nMovement'), 'fontsize', 15)
    xlabel(sprintf('Time, sec'), 'fontsize', 15)
end % END IF

% xTick = [(trialSize/2 - buffer)  + trialSize*([1:numTrials]-1)]./fs;
% xTick = [xTick, xTick + trialSize*numTrials./fs, xTick + 2*trialSize*numTrials./fs];
% % xTick = [xTick, xTick + trialSize*numTrials./fs];



% set(gca, 'ytick', [1:numChans] * 10 - 5)
% set(gca, 'yticklabel', [flipud(c2e(channels(:), 'pns'))])
% set(gca, 'xtick', xTick)
% set(gca, 'xticklabel', [1:numTrials])

set(gca, 'ytick', [])
set(gca, 'yticklabel', '')
set(gca, 'xtick', [0, 1, 2, 3])
set(gca, 'xticklabel', [0, 1, 2, 3])




end % END FUNTION

% EOF
