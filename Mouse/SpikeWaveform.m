%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	SpikeWaveform
%		Kevin O'Neill
%		Mario Capecchi Lab
%		20131103
%		v0.1
%		PI: Naveen Nagarajan
%
%	Inputs:
%		SpikeData: Data loaded from the expYYYY-MM-DD_hh-mm-ss_SE.mat file
%											(nse2mat output)
%
%       saveData: Boolean. If true, the figures will be saved.
%
%       channel: Scalar. Selects the channel to be analyzed. If needed can
%                be converted into a vector by adding a simple FOR loop.
%
%       unit: scalar. Selects the spike unit that is desired to be plotted.
%                     Requires Spike Sorting through a third party program.
%
%       ADDITIONAL NOTE: Instead of using patches it may be desired to just
%       overlay multiple waveforms. In which case there needs to be a
%       boolean flag to handle the option and use of the command: 
%               plot(SpikeData.samples, 'k')
%
%	Outputs:
%		There are no function outputs. But there is a saved file output.
%		Figures
%
%	To Use:
%		Run function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = SpikeWaveform( SpikeData, saveData, channel, unit)

if nargin < 4
    fprintf('Not enough input arguments to SpikeWaveform. Quitting program. DEBUG_nargin\n');
    return
end %END IF

% 	unitIdx = SpikeData.units == unit;
unitIdx = ones(1, length(SpikeData.timeStamps{channel}));

meanWave  = mean(SpikeData.samples{channel}(:,:), 2);
errorWave = std(SpikeData.samples{channel}(:,:), 0, 2);
errorWave = errorWave/sqrt(numel(SpikeData.samples)); % Standard error
errorWave = 2*errorWave; % 2*standard error.

x = 1:32;
xx = [x, fliplr(x)];

data1 =  [[meanWave + errorWave]', fliplr([meanWave - errorWave]')];

data1 = data1 * SpikeData.ADBitVolt;
meanData = meanWave * SpikeData.ADBitVolt;

figure(channel)
hold on
lData = plot(x, meanData, 'b');
pData = patch(xx, data1, 1);
lData = plot(x, meanData, 'b', 'linewidth', 2);

% Plots scale bar. Need to make this a little more elegant 
scaleData = plot([20,20], [max(data1), max(data1)-2e-5], 'k', 'linewidth', 2);
hold off

set(pData, 'FaceColor', 'k')
set(pData, 'EdgeColor', 'none')
set(pData, 'FaceAlpha', 0.25)

xlim([1, 32])
ylim([1.1*min(data1), 1.1*max(data1)])
%***** Hard Code ylim values if needed *****%

title('')
xlabel('')
ylabel('')

set(gca, 'XTick', [0])
set(gca, 'XTickLabel',{''})

set(gca, 'YTick', [0])
set(gca, 'YTickLabel',{''})

set(gca, 'Visible', 'Off')

text(21, max(data1)-1e-5, '2 \muV', 'FontWeight', 'bold')

set(gcf, 'Units', 'inches')
set(gcf, 'Position',[2 2 2 2])
set(gcf, 'PaperUnits','inches','PaperPosition',[2 2 2 2])

shortName = ['WaveForm_', SpikeData.expDate, '_chan', num2str(channel), '_unit', num2str(0)];

% if saveData
%     print('-dpng', ['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.png'], '-r100');
%     saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.fig']);
%     %		saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.png']);
%     saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.eps'],'epsc2');
% end % END IF

end % END FUNCTION