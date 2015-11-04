function [ output_args ] = RasterPlot( spikeTimes, varargin )

numIn = nargin;

%%%%% INIT %%%%%
fs         = 10000; % default sampling rate
flagLabel  = 1;     % boolean to show labels
spikeWidth = 1;     % spike thickness
spikeColor = 'k';   % spike color


switch numIn
    case 1 % Use defaults on Spike Times
        figure;
    case 2 % Alters sampling rate
        figure;
        fs = varargin{1};
    case 3 % Alters labels
        figure;
        fs        = varargin{1};
        flagLabel = varargin{2};
    case 4 % Spike Width
        figure;
        fs         = varargin{1};
        flagLabel  = varargin{2};
        spikeWidth = varargin{3};
    case 5 % Spike Color
        figure;
        fs         = varargin{1};
        flagLabel  = varargin{2};
        spikeWidth = varargin{3};
        spikeColor = varargin{4}; 
    otherwise
        error ('Invalid Arguments');
end

%%%%% MAIN %%%%%

spikeNum = length(spikeTimes);
xx = ones(1, 3 * spikeNum) * nan;
yy = ones(1, 3 * spikeNum) * nan;

xx(2:3:3*spikeNum) = spikeTimes .* (1000/fs);
xx(1:3:3*spikeNum) = xx(2:3:3*spikeNum) - eps(1);
xx(3:3:3*spikeNum) = xx(2:3:3*spikeNum) + eps(1);

yy(2:3:3*spikeNum) = 1;
yy(1:3:3*spikeNum) = 0;

plot(xx, yy, spikeColor, 'linewidth', spikeWidth)

if flagLabel
    ylabel('Channel')
    xlabel('Time, ms')
end % END IF

set(gca, 'ytick', [])

end % END FUNCTION
% EOF
