function [] = RasterComp( trialStruct, trialTypes, trials, channels, trialLen, spikeTimes, varargin )

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
    case 6 % Use defaults on Spike Times

    case 7 % Alters sampling rate
        figure;
        fs = varargin{1};
    case 8 % Alters labels
        figure;
        fs        = varargin{1};
        flagLabel = varargin{2};
    case 9 % Spike Width
        figure;
        fs         = varargin{1};
        flagLabel  = varargin{2};
        spikeWidth = varargin{3};
    case 10 % Spike Color
        figure;
        fs         = varargin{1};
        flagLabel  = varargin{2};
        spikeWidth = varargin{3};
        spikeColor = varargin{4};
    otherwise
        error ('Invalid Arguments');
end

%%%%% MAIN %%%%%

% Parses information

numTrials    = length(trials);
numChans     = size(channels, 1);
numMoves     = length(trialTypes);


% Loops over the number of channels
for i = 1:numChans
    
    H = figure(i);
    for j = 1:numMoves
        
        trialType    = trialTypes{j};
        rasterTrials = eval(['trialStruct.', trialType,'(trials,:)']);
        
        % Loops over the number of trials
        for k = 1:numTrials
            
%             if j==3 && k==20 && i==17
%                 disp([k, i, j])
%             end

            % Size of each trial to plot (Samples)
            trialSize = trialLen;
            startTime = rasterTrials(k,4);
            
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
                
                yy(2:3:3*spikeNum) = (numMoves - (j-1)) * 10 - trialGap;
                yy(1:3:3*spikeNum) = (numMoves - (j-1)-1) * 10;
                hold on
                plot((xx), yy, spikeColor{j})
                hold off
            end % END IF
        end % END FOR
    end % END FOR
    
    % Styles plot
    xlim([0, 3])
    ylim([0, numMoves*10])
    
    if flagLabel
        ylabel('Movement Type')
        xlabel('Time, sec')
        title(sprintf('Electrode: %d', c2e(channels(i), 'pns')))
    end % END IF
    
    % xTick = [(trialSize/2 - buffer)  + trialSize*([1:numTrials]-1)]./fs;
    % xTick = [xTick, xTick + trialSize*numTrials./fs, xTick + 2*trialSize*numTrials./fs];
    % xTick = [xTick, xTick + trialSize*numTrials./fs];
    
    set(gca, 'ytick', [5:10:30])
    set(gca, 'yticklabel', repmat({''}, 1,3))
    set(gca, 'xtick', [0:3])
    set(gca, 'xticklabel', [0:3])
    
    text(-0.1, 25, 'Hello', 'HorizontalAlignment', 'Center', 'Rotation', 90)
    text(-0.1, 15, 'Every', 'HorizontalAlignment', 'Center', 'Rotation', 90)
    text(-0.1, 5, 'Body', 'HorizontalAlignment', 'Center', 'Rotation', 90)
      
end % END FOR

end % END FUNTION

% EOF
