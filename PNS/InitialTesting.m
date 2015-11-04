%% Init %%

load '20130415-120054-001.mat'

%% Waveform %%

waveform = NEV.Data.Spikes.Waveform;
waveformAvg = mean(waveform, 2);

plot([1:length(waveformAvg)]./30, waveformAvg)
title('Averave Waveform')
xlabel('Time, ms')
ylabel('Voltage, mV')

%% SpikeTimes %%

% load '20130415-120054-001.mat'
rasterMat = zeros(96, 20);
fieldList = fieldnames(trialStruct);
spikeTimes = NEV.Data.Spikes.TimeStamp;
channelTimes = NEV.Data.Spikes.Electrode;

for i = 1:96
    spikeChannel{i} = spikeTimes(find(channelTimes == i));
end % END FOR
spikeChannel = spikeChannel';

for i = 1:size(fieldList, 2)
    rasterMat = zeros(96, 20);
    for k = 1:10
        figure (k)
        for j = 1:10
            if j+10*(k-1) < 97
                H = subplot(2,5,j);
                rasterMat(j+10*(k-1), :) = RasterPlotTrial(trialStruct, fieldList{i}, [1:20],  60000, spikeChannel{j+10*(k-1)}, H);
                title(['Channel: ', num2str(j+10*(k-1))])
            end % END IF
        end % END FOR
    end % END FOR
    rasterMat = [[1:96]', sum(rasterMat,2)];
    rasterMat = flipud(sortrows(rasterMat,2));
    eval(['rasterStruct.', fieldList{i},' = rasterMat(1:10, :);'])
end % END FOR



rasterStruct.T00000200 = rasterMat(1:10, :);
RasterPlotTrial(trialStruct, 'T00020000', [1:20],  60000, spikeChannel{2})
% RasterPlotTrial(trialStruct, 'T00020000', [1:40],  45000, spikeTimes)


%T00000000 : 54, 49, 43, 23, 21, 8
%T00000002 : 96, 54, 49, 43, 29, 23, 21, 8
%T00000020 : 96, 63, 54, 49, 43, 29, 23, 21, 8
%T00000200 : 96, 63,54,49, 43, 48, 23, 21, 8
%T00001000 : 72, 54, 49, 48, 43, 21, 8, 2
%T00020000 : 96,63, 54, 49, 43, 29, 23, 21, 19, 8, 2


% RasterPlot(spikeChannel{3}(1:300), 30000)
% 
% boxWin = BoxcarWindow(100, 30000);
% plot([1:length(boxWin)]./30, boxWin)

% yTest = zeros(1,spikeChannel{3}(25)+1);
% yTest(spikeChannel{3}(1:25)) = 1;
% 
% yTestConv = conv(yTest, boxWin);


% % % Plotting
% % figure (2)
% % plot([1:length(yTestConv)]./30, yTestConv)
% % 
% % hold on
% % 
% % windowSize = 50 * 30;
% % filtY = filter(ones(1,windowSize)/windowSize,1,yTestConv);
% % plot([1:length(filtY)]./30, filtY, 'r')
% % 
% % hold off
% % 
% % title('Action Potential Boxcar and Firing Rate')
% % xlabel('Time, ms')
% % ylabel('Firing Rate (APs/bin)')
% % legend('Boxcar, 100ms bin', 'Firing Rate, 50ms window')

%% Markers %%

% humanData = readNSxMat(...
%     ['D:\Kevin\PNS\P201301\20130415-120054\', fileName]);

load '20130415-120054-001.mat'
paramStruct = parseNEV_FingerPress_Kevin(NEV);
trialStruct = parseParamStruct(paramStruct);


%% Fatigue %%

clc, clear

load '20130415-120054-001_Sorted_TD.mat'
paramStruct = parseNEV_FingerPress_Kevin(NEV);
trialStruct = parseParamStruct(paramStruct);

%% 

spikeTimes = NEV.Data.Spikes.TimeStamp;
channelTimes = NEV.Data.Spikes.Electrode;
fieldList = fieldnames(trialStruct);

for i = 1:96
    spikeChannel{i} = spikeTimes(find(channelTimes == i));
end % END FOR
spikeChannel = spikeChannel';

H = subplot(1,1,1);

for k = 1:length(fieldList)
    H = subplot(1,1,1);
    RasterPlotElectrode(trialStruct, fieldList{k}, [1:20],  corrStrut, 90000, spikeChannel, H);
    title(sprintf('Task: %s\nChannels 1-96 for 20 trials',fieldList{k}))
    
    saveas(gcf,['D:\Kevin\PNS\P201301\Figures\RasterTime_', fieldList{k}, '.fig']);
    saveas(gcf,['D:\Kevin\PNS\P201301\Figures\RasterTime_', fieldList{k}, '.png']);
    saveas(gcf,['D:\Kevin\PNS\P201301\Figures\RasterTime_', fieldList{k}, '.eps'],'epsc2');
    
    close(gcf)

end % END FOR

% EOF

%% Correletion

load '20130415-120054-001.mat'
paramStruct = parseNEV_FingerPress_Kevin(NEV);
trialStruct = parseParamStruct(paramStruct);

%

% corrOut = PNSCorr( 'D:\Kevin\PNS\P201301\20130415-120054\20130415-120054-001', trialStruct, 60);
corrOut = PNSCorr( 'C:\CodeRepo\Lab\PNS\P201301\20130415-120054\20130415-120054-001', trialStruct, 60);

corrOut(:,15) = 1:96;

fieldList = fieldnames(trialStruct);

for k = 1:size(fieldList, 1)
    corrOut = flipud(sortrows(corrOut, k));
    eval(['corrStruct.',fieldList{k},'= corrOut([1:5], [k, end]);'])
end % END FOR


% HeatMap(corrOut, 'Symmetric', false, 'Colormap', jet)
corrOut2 = corrOut(:,1:end-1);
heatMax = max(corrOut2(:));

imgHeat = zeros(size(corrOut2, 1), size(corrOut2,2), 3);
mapColor = jet;

for i = 1:size(corrOut2, 1)
    for j = 1:size(corrOut2,2)
    
    imgHeat(i,j,:) = mapColor(ceil(corrOut2(i,j)./heatMax *  64), :);
    
    end % END FOR
    
end % END FOR

imagesc(imgHeat)

spikeTimes = NEV.Data.Spikes.TimeStamp;
channelTimes = NEV.Data.Spikes.Electrode;
fieldList = fieldnames(trialStruct);

for i = 1:96
    spikeChannel{i} = spikeTimes(find(channelTimes == i));
end % END FOR
spikeChannel = spikeChannel';

% H = subplot(1,1,1);

%%
for k = 1:length(fieldList)
    H = subplot(1,1,1);
    RasterPlotCorr(trialStruct, fieldList{k}, [1:20],  channels, 90000, spikeChannel, H);
    title(sprintf('Task: %s\nHighly Correleated Channels for 20 trials',fieldList{k}))
    
    saveas(gcf,['D:\Kevin\PNS\P201301\Figures\RasterCorr_', fieldList{k}, '.fig']);
    saveas(gcf,['D:\Kevin\PNS\P201301\Figures\RasterCorr_', fieldList{k}, '.png']);
    saveas(gcf,['D:\Kevin\PNS\P201301\Figures\RasterCorr_', fieldList{k}, '.eps'],'epsc2');
    
    close(gcf)

end % END FOR


%% Correllation Comparison %%

load '20130415-120054-001.mat'
paramStruct = parseNEV_FingerPress_Kevin(NEV);
trialStruct = parseParamStruct(paramStruct);
electrodes = [79, 6, 97, 36, 89, 26, 88, 70, 13, 44, 34, 84, 41, 14, 42, 12, 45]';
channels = e2c(electrodes, 'PNS');

% channels1 = elec2chan(electrodes);
channels1 = e2c(electrodes);
channels2 = e2c(electrodes, 'pns');


%%

spikeTimes = NEV.Data.Spikes.TimeStamp;
channelTimes = NEV.Data.Spikes.Electrode;
fieldList = fieldnames(trialStruct);

% fieldList2 = {'T00000200', 'T00000020', 'T00000002'}; % Intrinsic
fieldList2 = {'T00200000', 'T00020000', 'T00002000'}; % M,R,L Flexion
% fieldList2 = {'T20000000', 'T02000000', 'T00200000'}; % T,I,M Flexion

for i = 1:96
    spikeChannel{i,1} = spikeTimes(find(channelTimes == i));
end % END FOR

%%
for k = 1:12
    H = subplot(1,1,1);
    RasterPlotCorrComparison(trialStruct, fieldList(k:k+2), [1:20],  channels2, 90000, spikeChannel, H);
%     RasterPlotCorrComparison(trialStruct, fieldList2, [1:20],  channels2, 90000, spikeChannel, H);
    title(sprintf('Tasks: %s, %s, %s\n20 Trials for Electrodes: [79, 6, 97, 36, 89, 26, 88, 70, 13, 44, 34, 84, 41, 14, 42, 12, 45]',...
        fieldList{k}, fieldList{k+1}, fieldList{k+2}));

    set(gcf, 'Units', 'inches')
    set(gcf, 'Position', [1 1 16 6])
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 6])
    
%     %%%% No labels %%%%
%     print('-dpng', ['D:\Kevin\PNS\P201301\Figures\RasterComp_', num2str(k), '.png'], '-r100');
%     saveas(gcf,['D:\Kevin\PNS\P201301\Figures\RasterComp_', num2str(k), '.fig']);
% %     saveas(gcf,['D:\Kevin\PNS\P201301\Figures\RasterComp_', num2str(k), '.png']);
%     saveas(gcf,['D:\Kevin\PNS\P201301\Figures\RasterComp_', num2str(k), '.eps'],'epsc2');

%     close(gcf)    
end

%% More Plotting

    H = subplot(1,1,1);
%     RasterPlotCorrComparison(trialStruct, fieldList(k:k+2), [1:20],  channels2, 90000, spikeChannel, H);
    RasterPlotCorrComparison2(trialStruct, fieldList2, [1:20],  channels2, 90000, spikeChannel, H);
    title('Electrode: 26', 'FontSize', 15)
    
        set(gcf, 'Units', 'inches')
    set(gcf, 'Position',[0 0 6 2])
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 2])
    
     set(gca,'FontSize',15)
    
    %     title(sprintf('Tasks: %s, %s, %s\n20 Trials for Electrodes: [79, 6, 97, 36, 89, 26, 88, 70, 13, 44, 34, 84, 41, 14, 42, 12, 45]',...
%         fieldList{k}, fieldList{k+1}, fieldList{k+2}));

%     title(sprintf('Tasks: %s, %s, %s\n20 Trials for Electrodes: [79, 6, 97, 36, 89, 26, 88, 70, 13, 44, 34, 84, 41, 14, 42, 12, 45]',...
%         'Index Intrinsic', 'Ring Intrinsic', 'Little Intrinsic'));

%     title(sprintf('Tasks: %s, %s, %s\n20 Trials for Electrodes: [79, 6, 97, 36, 89, 26, 88, 70, 13, 44, 34, 84, 41, 14, 42, 12, 45]',...
%         'Middle Flexion', 'Ring Flexion', 'Little Flexion'));

% legend('Middle Flex', 'Ring Flex', 'Little Flex')
%%
title(sprintf('Middle Flex                                      Ring Flex                                      Little Flex'), 'fontsize', 15)

    set(gcf, 'Units', 'inches')
    set(gcf, 'Position', [1 1 12 6])
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 6])
    
    
    %%%% No labels %%%%
%     print('-dpng', ['D:\Kevin\PNS\P201301\Figures\RasterComp_', num2str(k), '.png'], '-r100');
%     saveas(gcf,['D:\Kevin\PNS\P201301\Figures\RasterComp_', num2str(k), '.fig']);
% %     saveas(gcf,['D:\Kevin\PNS\P201301\Figures\RasterComp_', num2str(k), '.png']);
%     saveas(gcf,['D:\Kevin\PNS\P201301\Figures\RasterComp_', num2str(k), '.eps'],'epsc2');

    %%%% Labels (eg: Index Flex....) %%%
%     print('-dpng', ['D:\Kevin\PNS\P201301\Figures\RasterComp_Flexion_MRL.png'], '-r100');
%     saveas(gcf,['D:\Kevin\PNS\P201301\Figures\RasterComp_Flexion_MRL.fig']);
% %     saveas(gcf,['D:\Kevin\PNS\P201301\Figures\RasterComp_Intrinsic.png']);
%     saveas(gcf,['D:\Kevin\PNS\P201301\Figures\RasterComp_Flexion_MRL.eps'],'epsc2');
    
%     close(gcf)

%% Firing Rate Comparison

load '20130415-120054-001.mat'
paramStruct = parseNEV_FingerPress_Kevin(NEV);
trialStruct = parseParamStruct(paramStruct);
electrodes = [79, 6, 97, 36, 89, 26, 88, 70, 13, 44, 34, 84, 41, 14, 42, 12, 45]';
channels = e2c(electrodes, 'pns');

spikeTimes = NEV.Data.Spikes.TimeStamp;
channelTimes = NEV.Data.Spikes.Electrode;
fieldList = fieldnames(trialStruct);

for i = 1:96
    spikeChannel{i,1} = spikeTimes(find(channelTimes == i));
end % END FOR

fieldList2 = {'T00200000', 'T00020000', 'T00002000'}; % M,R,L Flexion
fieldWord  = {'Middle Flex', 'Ring Flex', 'Little Flex'}; % M,R,L Flexion

% fieldList2 = {'T20000000', 'T02000000', 'T00200000'}; % T,I,M Flexion
% fieldWord  = {'Thumb Flex', 'Index Flex', 'Middle Flex'}; % T,I,M Flexion

% fieldList2 = {'T00000200', 'T00000020', 'T00000002'}; % I,R,L Intrinsic
% fieldWord  = {'Index Intrinsic', 'Ring Intrinsic', 'Little Intrinsic'}; % I,R,L Intrinsic


%%

for k = 1:length(fieldList2)%2:4
    subplot(1,3,k)
    firingRate = FiringRate(trialStruct, fieldList2{k}, [1:20],  channels, 90000, spikeChannel);
    

    boxplot(firingRate)
%     title(sprintf('Task: %s\nBox and Whisker plot of AP firing rate over 20 trials', fieldList2{k}))
    title(sprintf('Task: %s\nBox and Whisker plot of AP firing rate over 20 trials', fieldWord{k}))

    xlabel('Electrodes')
    ylabel('Action Potential Firing Rate, Hz')
    ylim([0, 25])
    set(gca, 'xtick', 1:17)
    set(gca, 'xticklabel', [c2e(channels(:)', 'pns')])
       
end % END FOR

set(gcf, 'Units', 'inches')
set(gcf, 'Position', [1 1 16 6])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 6])

%%% No labels %%%%
% print('-dpng', ['D:\Kevin\PNS\P201301\Figures\FiringRate_IRL_Intrinsic.png'], '-r100');
% saveas(gcf,['D:\Kevin\PNS\P201301\Figures\FiringRate_IRL_Intrinsic.fig']);
% %     saveas(gcf,['D:\Kevin\PNS\P201301\Figures\FiringRate_IRL_Intrinsic.png']);
% saveas(gcf,['D:\Kevin\PNS\P201301\Figures\FiringRate_IRL_Intrinsic.eps'],'epsc2');


%%  Firing Rate Continuous

load '20130415-120054-001.mat'
paramStruct = parseNEV_FingerPress_Kevin(NEV);
trialStruct = parseParamStruct(paramStruct);
electrodes = [79, 6, 97, 36, 89, 26, 88, 70, 13, 44, 34, 84, 41, 14, 42, 12, 45]';
channels = e2c(electrodes, 'pns');

spikeTimes = NEV.Data.Spikes.TimeStamp;
channelTimes = NEV.Data.Spikes.Electrode;
fieldList = fieldnames(trialStruct);

for i = 1:96
    spikeChannel{i,1} = spikeTimes(find(channelTimes == i));
end % END FOR

% fieldList2 = {'T00200000', 'T00020000', 'T00002000'}; % M,R,L Flexion
% fieldWord  = {'Middle Flex', 'Ring Flex', 'Little Flex'}; % M,R,L Flexion
% compList = 'MRLflex';

fieldList2 = {'T20000000', 'T02000000', 'T00200000'}; % T,I,M Flexion
fieldWord  = {'Thumb Flex', 'Index Flex', 'Middle Flex'}; % T,I,M Flexion
compList = 'TIMflex';

% fieldList2 = {'T00000200', 'T00000020', 'T00000002'}; % I,R,L Intrinsic
% fieldWord  = {'Index Intrinsic', 'Ring Intrinsic', 'Little Intrinsic'}; % I,R,L Intrinsic
% compList = 'IRLabd';

%%
rateMat = FiringRateCont(trialStruct, fieldList2{1}, [1:20],  channels, 90000, spikeChannel);

rateMat2 = mean(rateMat, 3);
x=linspace(-0.5, 3.5, 120000);

%%
for i = 1:3
    subplot(3, 1, i)
    xlabel('Time, sec')
    ylabel('AP Firing Rate, Hz')
    
    rateMat = FiringRateCont(trialStruct, fieldList2{i}, [1:20],  channels, 90000, spikeChannel);
    
    rateMat2 = mean(rateMat, 3);
    x=linspace(-0.5, 3.5, 120000);
    
    for k = 1:17

        subplot(3, 17, k + (i-1)*17)        
        
        plot(x, rateMat2(k,:).*30000)
        if i == 1
            title(['Elec: ', num2str(electrodes(k))])
        else
            title('')
        end % END IF
        
        xlim([-0.5, 3.5])
        %     ylim([0, 8*10^-4])
        ylim([0, 30])
        
        if k == 1
            set(gca, 'ytick', [0, 15, 30])
            set(gca, 'yticklabel', [0, 15, 30])
            set(gca, 'xtick', [0, 1, 2, 3])
            set(gca, 'xticklabel', [0, 1, 2, 3])
            
            p1 = get(gca, 'position');
        else
            set(gca, 'ytick', [0, 15, 30])
            set(gca, 'yticklabel', ['', '', ''])
            set(gca, 'xtick', [0, 1, 2, 3])
            set(gca, 'xticklabel', [0, 1, 2, 3])
        end % END FOR
        
        if k == 17 && i == 1
            p2 = get(gca, 'position');
            width = p2(1) + p2(3) - p1(1);
            h3 = axes('position',[p1(1), p1(2), width, p1(4)],'visible','off');
            hLabelY = ylabel(sprintf('Trial: %s\nAP Firing Rate, Hz', fieldWord{i}),'visible','on');
            hLabelX = xlabel('Time, sec','visible','on');
            
            h3 = axes('position',[p1(1), p1(2), width, p1(4)+0.04],'visible','off');
            hTitle = title('Average Continuous AP Firing Rate Across 20 Trials','visible','on');
        elseif k == 17
            p2 = get(gca, 'position');
            width = p2(1) + p2(3) - p1(1);
            h3 = axes('position',[p1(1), p1(2), width, p1(4)],'visible','off');
            hLabelY = ylabel(sprintf('Trial: %s\nAP Firing Rate, Hz', fieldWord{i}),'visible','on');
            hLabelX = xlabel('Time, sec','visible','on');
        end % END IF
    end % END FOR
end % END FOR

set(gcf, 'Units', 'inches')
set(gcf, 'Position', [1 1 16 6])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 6])

%%% No labels %%%%
print('-dpng', ['D:\Kevin\PNS\P201301\Figures\FiringRateCont_', compList, '.png'], '-r100');
saveas(gcf,['D:\Kevin\PNS\P201301\Figures\FiringRateCont_', compList, '.fig']);
%     saveas(gcf,['D:\Kevin\PNS\P201301\Figures\FiringRateCont_', compList, '.png']);
saveas(gcf,['D:\Kevin\PNS\P201301\Figures\FiringRateCont_', compList, '.eps'],'epsc2');


%%

for k = 1:length(fieldList2)%2:4
    subplot(1,3,k)
    firingRate = FiringRate(trialStruct, fieldList2{k}, [1:20],  channels, 90000, spikeChannel);
    

    boxplot(firingRate)
%     title(sprintf('Task: %s\nBox and Whisker plot of AP firing rate over 20 trials', fieldList2{k}))
    title(sprintf('Task: %s\nBox and Whisker plot of AP firing rate over 20 trials', fieldWord{k}))

    xlabel('Electrodes')
    ylabel('Action Potential Firing Rate, Hz')
    ylim([0, 25])
    set(gca, 'xtick', 1:17)
    set(gca, 'xticklabel', [c2e(channels(:)', 'pns')])
       
end % END FOR

set(gcf, 'Units', 'inches')
set(gcf, 'Position', [1 1 16 6])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 6])

%%% No labels %%%%
print('-dpng', ['D:\Kevin\PNS\P201301\Figures\FiringRate_IRL_Intrinsic.png'], '-r100');
saveas(gcf,['D:\Kevin\PNS\P201301\Figures\FiringRate_IRL_Intrinsic.fig']);
%     saveas(gcf,['D:\Kevin\PNS\P201301\Figures\FiringRate_IRL_Intrinsic.png']);
saveas(gcf,['D:\Kevin\PNS\P201301\Figures\FiringRate_IRL_Intrinsic.eps'],'epsc2');

%% PSTH

load '20130415-120054-001_Sorted_TD.mat'
% load '20130415-120054-001.mat'
paramStruct = parseNEV_FingerPress_Kevin(NEV);
trialStruct = parseParamStruct(paramStruct);
electrodes = [79, 6, 97, 36, 89, 26, 88, 70, 13, 44, 34, 84, 41, 14, 42, 12, 45]';
channels = e2c(electrodes, 'pns');

spikeTimes = NEV.Data.Spikes.TimeStamp(find(NEV.Data.Spikes.Unit == 1));
channelTimes = NEV.Data.Spikes.Electrode(find(NEV.Data.Spikes.Unit == 1));
fieldList = fieldnames(trialStruct);

for i = 1:96
    spikeChannel{i,1} = spikeTimes(find(channelTimes == i));
end % END FOR

fieldList2 = {'T00200000', 'T00020000', 'T00002000'}; % M,R,L Flexion
fieldWord  = {'Middle Flex', 'Ring Flex', 'Little Flex'}; % M,R,L Flexion
compList = 'MRLflex';

% fieldList2 = {'T20000000', 'T02000000', 'T00200000'}; % T,I,M Flexion
% fieldWord  = {'Thumb Flex', 'Index Flex', 'Middle Flex'}; % T,I,M Flexion
% compList = 'TIMflex';

% fieldList2 = {'T00000200', 'T00000020', 'T00000002'}; % I,R,L Intrinsic
% fieldWord  = {'Index Intrinsic', 'Ring Intrinsic', 'Little Intrinsic'}; % I,R,L Intrinsic
% compList = 'IRLabd';


%%
figure (1)
rateMat = FiringRatePSTH(trialStruct, fieldList2{1}, [1:20],  channels, 90000, spikeChannel);

rateMat2 = rateMat;
rateMat2(:,:,:,2) = FiringRatePSTH(trialStruct, fieldList2{2}, [1:20],  channels, 90000, spikeChannel);
rateMat2(:,:,:,3) = FiringRatePSTH(trialStruct, fieldList2{3}, [1:20],  channels, 90000, spikeChannel);

% close(gcf)

x=linspace(-2, 5, length(rateMat2(:,1,1,1)));
xx = [x, fliplr(x)];

%%
% figure(2)
for k = 1:17
%     subplot(3, 6, k)
figure(k)
    dataTop1 = rateMat2(:,k,1,1)+rateMat2(:,k,2,1);
    dataBot1 = rateMat2(:,k,1,1)-rateMat2(:,k,2,1);
    
    dataTop2 = rateMat2(:,k,1,2)+rateMat2(:,k,2,2);
    dataBot2 = rateMat2(:,k,1,2)-rateMat2(:,k,2,2);
    
    dataTop3 = rateMat2(:,k,1,3)+rateMat2(:,k,2,3);
    dataBot3 = rateMat2(:,k,1,3)-rateMat2(:,k,2,3);
    
    approxArea = [sum(dataTop1) - sum(dataBot1); sum(dataTop2) - sum(dataBot2); sum(dataTop3) - sum(dataBot3)];
    approxArea = [approxArea, [1;2;3]];
    
    approxArea = flipud(sortrows(approxArea, 1));
    
    data1 = [[rateMat2(:,k,1,1)+rateMat2(:,k,2,1)]', fliplr([rateMat2(:,k,1,1)-rateMat2(:,k,2,1)]')];
    data2 = [[rateMat2(:,k,1,2)+rateMat2(:,k,2,2)]', fliplr([rateMat2(:,k,1,2)-rateMat2(:,k,2,2)]')];
    data3 = [[rateMat2(:,k,1,3)+rateMat2(:,k,2,3)]', fliplr([rateMat2(:,k,1,3)-rateMat2(:,k,2,3)]')];
    
    hold on
%     pTop1 = patch(x,rateMat2(:,k,1,1)+rateMat2(:,k,2,1), 1);
%     pTop2 = patch(x,rateMat2(:,k,1,2)+rateMat2(:,k,2,2), 1);
%     pTop3 = patch(x,rateMat2(:,k,1,3)+rateMat2(:,k,2,3), 1);

    tmpR = patch([-2, -2, -1, -1], [0, -1, -1, 0], 1);
    tmpG = patch([-2, -2, -1, -1], [0, -1, -1, 0], 1);
    tmpB = patch([-2, -2, -1, -1], [0, -1, -1, 0], 1);
    
    set(tmpR, 'FaceColor', 'r');
    set(tmpG, 'FaceColor', 'g');
    set(tmpB, 'FaceColor', 'b');

%     set(pTop1, 'FaceColor', 'r'); set(pTop2, 'FaceColor', 'g'); set(pTop3, 'FaceColor', 'b');
%     set(pTop1, 'EdgeColor', 'none'); set(pTop2, 'EdgeColor', 'none'); set(pTop3, 'EdgeColor', 'none');

    for i = 1:3
%         eval(sprintf('pTop%d = patch(x,dataTop%d, 1);', approxArea(i,2), approxArea(i,2)))
%         eval(sprintf('pBot%d = patch(x,dataBot%d, 1);', approxArea(i,2), approxArea(i,2)))
        eval(sprintf('pData%d = patch(xx,data%d, 1);', approxArea(i,2), approxArea(i,2)))
        eval(sprintf('lData%d = plot(x,rateMat2(:,k,1,%d));', approxArea(i,2), approxArea(i,2)))

    
    end % END FOR
%     pTop1 = patch(x,rateMat2(:,k,1,1)+rateMat2(:,k,2,1), 1);
%     pBot1 = patch(x,rateMat2(:,k,1,1)-rateMat2(:,k,2,1), 1); 
%     
%     pTop2 = patch(x,rateMat2(:,k,1,2)+rateMat2(:,k,2,2), 1);
%     pBot2 = patch(x,rateMat2(:,k,1,2)-rateMat2(:,k,2,2), 1); 
%     
%     pTop3 = patch(x,ra teMat2(:,k,1,3)+rateMat2(:,k,2,3), 1);
%     pBot3 = patch(x,rateMat2(:,k,1,3)-rateMat2(:,k,2,3), 1);
    
%     set(pTop1, 'FaceColor', 'r');     set(pBot1, 'FaceColor', 'w');
%     set(pTop2, 'FaceColor', 'g');     set(pBot2, 'FaceColor', 'w');
%     set(pTop3, 'FaceColor', 'b');     set(pBot3, 'FaceColor', 'w');
    
%     set(pTop1, 'EdgeColor', 'none'); set(pTop2, 'EdgeColor', 'none'); set(pTop3, 'EdgeColor', 'none');
%     set(pBot1, 'EdgeColor', 'none'); set(pBot2, 'EdgeColor', 'none'); set(pBot3, 'EdgeColor', 'none');


    set(pData1, 'FaceColor', 'r');
    set(pData2, 'FaceColor', 'g');
    set(pData3, 'FaceColor', 'b');
    
    set(pData1, 'FaceAlpha', 0.25);
    set(pData2, 'FaceAlpha', 0.25);
    set(pData3, 'FaceAlpha', 0.25);
    
    set(pData1, 'EdgeColor', 'none');
    set(pData2, 'EdgeColor', 'none');
    set(pData3, 'EdgeColor', 'none');

    set(lData1, 'Color', 'r');
    set(lData2, 'Color', 'g');
    set(lData3, 'Color', 'b');
    
    set(lData1, 'linewidth', 2.75);
    set(lData2, 'linewidth', 2.75);
    set(lData3, 'linewidth', 2.75);
    
    plot(x, zeros(length(dataTop1)), 'k')
    hold off
    
%     title(sprintf('PSTH Comparison\nBetween Movement Types\nElectrode: %d', electrodes(k)))
    title(sprintf('Electrode: %d', electrodes(k)), 'fontsize', 15)

    ylabel('Firing Rate, Hz', 'fontsize', 15)
    xlabel('Time, sec', 'fontsize', 15)

        set(gcf, 'Units', 'inches')
    set(gcf, 'Position',[0 0 6 3])
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 3])

    xlim([0, 3])
    ylim([0, 40])
    
%           xlabel('Time, sec', 'font size', 15)
        p1 = get(gca, 'position');
        legend(fieldWord{:}, 'Location', [0.79, 0.8, 0.125, 0.0625])
        legend boxoff
        
        set(gca, 'ytick', [0, 10, 20, 30, 40])
        set(gca, 'yticklabel', [0, 10, 20, 30, 40])
        set(gca, 'xtick', [0, 1, 2, 3])
        set(gca, 'xticklabel', [0, 1, 2, 3])    
%  grid on
 
 set(gca,'FontSize',15)
 
 h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',15); 

h_ylabel = get(gca,'yLabel');
set(h_xlabel,'FontSize',15); 

    print('-dpng', ['C:\CodeRepo\Lab\PNS\P201301\Figures2\PSTH_V2_', num2str(k), '.png'], '-r100');
    saveas(gcf,['C:\CodeRepo\Lab\PNS\P201301\Figures2\PSTH_V2_', num2str(k), '.png']);

        
%     if k == 1
% %         set(gca, 'ytick', [0, 15, 30])
% %         set(gca, 'yticklabel', [0, 15, 30])
% %         set(gca, 'xtick', [0, 1, 2, 3])
% %         set(gca, 'xticklabel', [0, 1, 2, 3])
%         xlabel('Time, sec')
%         p1 = get(gca, 'position');
%         legend(fieldWord{:}, 'Location', [0.89, 0.8, 0.125, 0.0625])
%     else
% %         set(gca, 'ytick', [0, 15, 30])
% %         set(gca, 'yticklabel', ['', '', ''])
% %         set(gca, 'xtick', [0, 1, 2, 3])
% %         set(gca, 'xticklabel', [0, 1, 2, 3])
%     end % END FOR
%     
%     if k == 17 
%         p2 = get(gca, 'position');
%         width = p2(1) + p2(3) - p1(1);
%         h3 = axes('position',[p1(1), p1(2), width, p1(4)],'visible','off');
%         hLabelY = ylabel(sprintf('Action Potential Firing Rate, Hz'),'visible','on');
%         
%         h3 = axes('position',[p1(1), p1(2), width, p1(4)+0.04],'visible','off');
%         hTitle = title('PSTH Comparison Between Movement Types','visible','on');
%     end % END IF
end % END FOR
% 
%     set(gcf, 'Units', 'inches')
%     set(gcf, 'Position',[0 0 16 8])
%     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 8])
    
%     set(gcf, 'Units', 'inches')
%     set(gcf, 'Position',[0 0 6 4])
%     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 4])
% 
%     legend(fieldWord{:}, 'Location', [0.89, 0.8, 0.125, 0.0625])
% 
%     print('-dpng', ['D:\Kevin\PNS\P201301\Figures\PSTH_V2_', compList, '.png'], '-r100');
%     saveas(gcf,['D:\Kevin\PNS\P201301\Figures\PSTH_V2_', compList, '.fig']);
% % %     saveas(gcf,['D:\Kevin\PNS\P201301\Figures\PSTH_V2_', compList, '.png']);
%     saveas(gcf,['D:\Kevin\PNS\P201301\Figures\PSTH_V2_', compList, '.eps'],'epsc2');
%     
% %     close(gcf)

%%
for i = 1:3
    subplot(3, 1, i)
    xlabel('Time, sec')
    ylabel('AP Firing Rate, Hz')
    
    
%     rateMat2 = mean(rateMat, 3);
    x=linspace(0, 3, length(rateMat(:,1,1)));
    
    for k = 1:17

        subplot(3, 17, k + (i-1)*17)        
        hold on
        plot(x, rateMat2(:,k,1,i), 'b')
%         plot(x, rateMat(:,k,2), 'r')
        hold off
        if i == 1
            title(['Elec: ', num2str(electrodes(k))])
        else
            title('')
        end % END IF
        
        xlim([0, 3])
        %     ylim([0, 8*10^-4])
        ylim([0, 30])
        
        if k == 1
            set(gca, 'ytick', [0, 15, 30])
            set(gca, 'yticklabel', [0, 15, 30])
            set(gca, 'xtick', [0, 1, 2, 3])
            set(gca, 'xticklabel', [0, 1, 2, 3])
            
            p1 = get(gca, 'position');
        else
            set(gca, 'ytick', [0, 15, 30])
            set(gca, 'yticklabel', ['', '', ''])
            set(gca, 'xtick', [0, 1, 2, 3])
            set(gca, 'xticklabel', [0, 1, 2, 3])
        end % END FOR
        
        if k == 17 && i == 1
            p2 = get(gca, 'position');
            width = p2(1) + p2(3) - p1(1);
            h3 = axes('position',[p1(1), p1(2), width, p1(4)],'visible','off');
            hLabelY = ylabel(sprintf('Trial: %s\nAP Firing Rate, Hz', fieldWord{i}),'visible','on');
            hLabelX = xlabel('Time, sec','visible','on');
            
            h3 = axes('position',[p1(1), p1(2), width, p1(4)+0.04],'visible','off');
            hTitle = title('Average Continuous AP Firing Rate Across 20 Trials','visible','on');
        elseif k == 17
            p2 = get(gca, 'position');
            width = p2(1) + p2(3) - p1(1);
            h3 = axes('position',[p1(1), p1(2), width, p1(4)],'visible','off');
            hLabelY = ylabel(sprintf('Trial: %s\nAP Firing Rate, Hz', fieldWord{i}),'visible','on');
            hLabelX = xlabel('Time, sec','visible','on');
        end % END IF
    end % END FOR
end % END FOR

set(gcf, 'Units', 'inches')
set(gcf, 'Position', [1 1 16 6])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 6])

%% Tuning Curve

load '20130415-120054-001_Sorted_TD.mat'
% load '20130415-120054-001.mat'
paramStruct = parseNEV_FingerPress_Kevin(NEV);
trialStruct = parseParamStruct(paramStruct);
electrodes = [79, 6, 97, 36, 89, 26, 88, 70, 13, 44, 34, 84, 41, 14, 42, 12, 45]';
channels = e2c(electrodes, 'pns');

spikeTimes = NEV.Data.Spikes.TimeStamp(find(NEV.Data.Spikes.Unit == 1));
channelTimes = NEV.Data.Spikes.Electrode(find(NEV.Data.Spikes.Unit == 1));
fieldList = fieldnames(trialStruct);

for i = 1:96
    spikeChannel{i,1} = spikeTimes(find(channelTimes == i));
end % END FOR

moveWord = {'Thumb Flex', 'Index Flex', 'Middle Flex', 'Ring Flex', 'Little Flex', ...
            'Index Abd.', 'Ring Abd.', 'Little Abd.', ...
            'Thumb Ext.', 'Index Ext.', 'Middle Ext.', 'Ring Ext.', 'Little Ext.', ...
            'No movement'};

fieldList2 = {'T00200000', 'T00020000', 'T00002000'}; % M,R,L Flexion
fieldWord  = {'Middle Flex', 'Ring Flex', 'Little Flex'}; % M,R,L Flexion

% fieldList2 = {'T20000000', 'T02000000', 'T00200000'}; % T,I,M Flexion
% fieldWord  = {'Thumb Flex', 'Index Flex', 'Middle Flex'}; % T,I,M Flexion

% fieldList2 = {'T00000200', 'T00000020', 'T00000002'}; % I,R,L Intrinsic
% fieldWord  = {'Index Intrinsic', 'Ring Intrinsic', 'Little Intrinsic'}; % I,R,L Intrinsic

%%

tuneC = TuningCurve(trialStruct, [1:20],  channels, 90000, spikeChannel);

meanC = mean(tuneC, 2);
stdC = std(tuneC, 0, 2)/sqrt(20);

for k = 1:17
    figure(k)
    hold on
    errorbar(1:14, meanC(k, 1, :), stdC(k, 1, :), 'k.')
    tmp = meanC(k, 1, :);
    plot(1:14, tmp(:), 'r.')
    hold off
    
    title(sprintf('Tuning Curve for\nAction Potential Firing Rate vs Movement Type\nElectrode: %d', electrodes(k)))
%     xlabel('Movement Type') % Replaced with 'text' command below
    ylabel('Average AP Firing Rate, Hz')
    xlim([0, 15])
    ylim([0, 30])
    
    set(gca, 'ytick', [0:5:30])
    set(gca, 'yticklabel', [0:5:30])
    set(gca, 'xtick', [1:14])
    set(gca, 'xticklabel', repmat({''}, 1, 14))
    
    for i = 1:14
        text(i, -0.5, moveWord{i}, 'HorizontalAlignment', 'Right', 'Rotation', 45)
    end % END FOR
    
    text(7, -6, 'Movement Type', 'HorizontalAlignment', 'Center')
    
    set(gcf, 'Units', 'inches')
    set(gcf, 'Position', [1 1 7 5])
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 5])
    
    aPos = get(gca, 'Position');
    set(gca, 'Position', [aPos(1), aPos(2)+0.05, aPos(3), aPos(4)-0.05])
    
    set(gca,'XGrid','off','YGrid','on','ZGrid','off')
    
    print('-dpng', ['D:\Kevin\PNS\P201301\Figures\TuningCurve_Chan', num2str(channels(k)), '_E', num2str(electrodes(k)), '.png'], '-r100');
    saveas(gcf,['D:\Kevin\PNS\P201301\Figures\TuningCurve_Chan', num2str(channels(k)), '_E', num2str(electrodes(k)), '.fig']);
%     saveas(gcf,['D:\Kevin\PNS\P201301\Figures\PSTH_Chan', num2str(channels(k)), '_E', num2str(electrodes(k)), '.png']);
    saveas(gcf,['D:\Kevin\PNS\P201301\Figures\TuningCurve_Chan', num2str(channels(k)), '_E', num2str(electrodes(k)), '.eps'],'epsc2');
    
    close(gcf)
    
end % END FOR


%% Rasters for PSTH


load '20130415-120054-001_Sorted_TD.mat'
% load '20130415-120054-001.mat'
paramStruct = parseNEV_FingerPress_Kevin(NEV);
trialStruct = parseParamStruct(paramStruct);
electrodes = [79, 6, 97, 36, 89, 26, 88, 70, 13, 44, 34, 84, 41, 14, 42, 12, 45]';
channels = e2c(electrodes, 'pns');

spikeTimes = NEV.Data.Spikes.TimeStamp(find(NEV.Data.Spikes.Unit == 1));
channelTimes = NEV.Data.Spikes.Electrode(find(NEV.Data.Spikes.Unit == 1));
fieldList = fieldnames(trialStruct);

for i = 1:96
    spikeChannel{i,1} = spikeTimes(find(channelTimes == i));
end % END FOR

fieldList2 = {'T00200000', 'T00020000', 'T00002000'}; % M,R,L Flexion
fieldWord  = {'Middle Flex', 'Ring Flex', 'Little Flex'}; % M,R,L Flexion
compList = 'MRLflex';

% fieldList2 = {'T20000000', 'T02000000', 'T00200000'}; % T,I,M Flexion
% fieldWord  = {'Thumb Flex', 'Index Flex', 'Middle Flex'}; % T,I,M Flexion
% compList = 'TIMflex';

% fieldList2 = {'T00000200', 'T00000020', 'T00000002'}; % I,R,L Intrinsic
% fieldWord  = {'Index Intrinsic', 'Ring Intrinsic', 'Little Intrinsic'}; % I,R,L Intrinsic
% compList = 'IRLabd';

%%

RasterComp(trialStruct, fieldList2, [10],  channels, 90000, spikeChannel);



%% Voting Matrix

load '20130415-120054-001_Sorted_TD.mat'
% load '20130415-120054-001.mat'
paramStruct = parseNEV_FingerPress_Kevin(NEV);
trialStruct = parseParamStruct(paramStruct);
electrodes = [79, 6, 97, 36, 89, 26, 88, 70, 13, 44, 34, 84, 41, 14, 42, 12, 45]';
channels = e2c(electrodes, 'pns');

spikeTimes = NEV.Data.Spikes.TimeStamp(find(NEV.Data.Spikes.Unit == 1));
channelTimes = NEV.Data.Spikes.Electrode(find(NEV.Data.Spikes.Unit == 1));
fieldList = fieldnames(trialStruct);

for i = 1:96
    spikeChannel{i,1} = spikeTimes(find(channelTimes == i));
end % END FOR

moveWord = {'Thumb Flex', 'Index Flex', 'Middle Flex', 'Ring Flex', 'Little Flex', ...
            'Index Abd.', 'Ring Abd.', 'Little Abd.', ...
            'Thumb Ext.', 'Index Ext.', 'Middle Ext.', 'Ring Ext.', 'Little Ext.', ...
            'No movement'};

% fieldList2 = {'T00200000', 'T00020000', 'T00002000'}; % M,R,L Flexion
% fieldWord  = {'Middle Flex', 'Ring Flex', 'Little Flex'}; % M,R,L Flexion

% fieldList2 = {'T20000000', 'T02000000', 'T00200000'}; % T,I,M Flexion
% fieldWord  = {'Thumb Flex', 'Index Flex', 'Middle Flex'}; % T,I,M Flexion

% fieldList2 = {'T00000200', 'T00000020', 'T00000002'}; % I,R,L Intrinsic
% fieldWord  = {'Index Intrinsic', 'Ring Intrinsic', 'Little Intrinsic'}; % I,R,L Intrinsic

%%

tuneC = TuningCurve(trialStruct, [1:20],  channels, 90000, spikeChannel);

meanC = mean(tuneC, 2);
meanC = permute(meanC, [1,3,2]);

stdC = std(tuneC, 0, 2)/sqrt(20);
stdC = permute(stdC, [1, 3, 2]);

weightMat = meanC/max(sum(meanC, 1));

weightMat = weightMat - repmat(mean(weightMat,2), [1, 14]);

weightMat(weightMat <= 0) = NaN;

HeatMap(weightMat, 'ColorMap', 'copper', 'Symmetric', 'false');


%% PCA Gating

load '20130415-120054-001_Sorted_TD.mat'
% load '20130415-120054-001.mat'
paramStruct = parseNEV_FingerPress_Kevin(NEV);
trialStruct = parseParamStruct(paramStruct);
electrodes = [79, 6, 97, 36, 89, 26, 88, 70, 13, 44, 34, 84, 41, 14, 42, 12, 45]';
channels = e2c(electrodes, 'pns');

spikeTimes = NEV.Data.Spikes.TimeStamp(find(NEV.Data.Spikes.Unit == 1));
channelTimes = NEV.Data.Spikes.Electrode(find(NEV.Data.Spikes.Unit == 1));
fieldList = fieldnames(trialStruct);

for i = 1:96
    spikeChannel{i,1} = spikeTimes(find(channelTimes == i));
end % END FOR

fieldList2 = {'T00200000', 'T00020000', 'T00002000'}; % M,R,L Flexion
fieldWord  = {'Middle Flex', 'Ring Flex', 'Little Flex'}; % M,R,L Flexion
compList = 'MRLflex';

% fieldList2 = {'T20000000', 'T02000000', 'T00200000'}; % T,I,M Flexion
% fieldWord  = {'Thumb Flex', 'Index Flex', 'Middle Flex'}; % T,I,M Flexion
% compList = 'TIMflex';

% fieldList2 = {'T00000200', 'T00000020', 'T00000002'}; % I,R,L Intrinsic
% fieldWord  = {'Index Intrinsic', 'Ring Intrinsic', 'Little Intrinsic'}; % I,R,L Intrinsic
% compList = 'IRLabd';

%%

PCA = calcPCA(6, 1:14);

PCA.Channels = e2c(PCA.Electrodes,'pns');

[V,D] = eig(cov(PCA.Results.TrainMat));
V = fliplr(V); D = fliplr(D); %organizing vectors so most variance is in left column
Var90 = find(cumsum(sum(D))/sum(D(:))<=0.9,1,'last'); %finding first set of vectors with 90 percent of the variance
if isempty(Var90)
    Var90 = 1;
end

% 75672785 80489300

%%
currTime = 75672785 - 2*30000;
L = zeros(1, 3000);
L2 = zeros(1,3000);
for k = 1:400
    L(k) = ClassifyPCA(trialStruct, spikeChannel, PCA, 2, currTime + 15000*k);
    disp(k);
    
    if k>=15
        LHist = [LHist(2:end), L(k)];
    else
        LHist = L(1:k);
    end
    
    L2(k) = mean(LHist);
end % END FOR

LTest = L;


trialTime = ceil((trialStruct.T10000000(2:23, 4) - (75672785 - 2*30000))./2000);
if trialTime(1) == 0
    trialTime(1) = 1;
end % END IF
LTrue = ones(1, 400) * 14;
for i = 1:length(trialTime)
    LTrue(trialTime(i):trialTime(i)+3*15) = 1;
end % END FOR

hold on
plot(LTrue, 'k')
plot(LTest, 'b')
hold off


% FTest = reshape(permute(Features(FeatureIdxs,e2c(Elects,'pns'),11:20,:),[3,1,4,2]),10*length(FeatureIdxs),[]);
% FTestID = reshape(permute(FeatureID(FeatureIdxs,e2c(Elects,'pns'),11:20,:),[3,1,4,2]),10*length(FeatureIdxs),[]);
% FTestID = FTestID(~all(isnan(FTest),2),:); FTest = FTest(~all(isnan(FTest),2),:);
% 
% 
% 
% TrainMat = [FTrain;RTrain];
% 
% 
% [V,D] = eig(cov(PCA.Results.TrainMat));
% V = fliplr(V); D = fliplr(D); %organizing vectors so most variance is in left column
% Var90 = find(cumsum(sum(D))/sum(D(:))<=0.9,1,'last'); %finding first set of vectors with 90 percent of the variance
% if isempty(Var90)
%     Var90 = 1;
% end
% 
% PCTrainMat = TrainMat*V(:,1:Var90); %projecting to pc space
% PCTestMat = TestMat*V(:,1:Var90);
% 
% L = classify(PCTestMat,PCTrainMat,TrainMatID(:,1)); %classifying

%% PCA Plotting

x = [0.5, 1.5, 0,...
    1.5, 2.5, 0,...
    2.5, 3.5, 0,...
    3.5, 4.5, 0,...
    4.5, 5.5, 0,...
    5.5, 6.5, 0,...
    6.5, 7.5, 0,...
    7.5, 8.5, 0,...
    8.5, 9.5, 0,...
    9.5, 10.5,0,...
    10.5, 11.5, 0,...
    11.5, 12.5, 0,...
    12.5, 13.5];
y = [1/1, 1/1, nan,...
    1/2, 1/2, nan,...
    1/3, 1/3, nan,...
    1/4, 1/4, nan,...
    1/5, 1/5, nan,...
    1/6, 1/6, nan,...
    1/7, 1/7, nan,...
    1/8, 1/8, nan,...
    1/9, 1/9, nan,...
    1/10, 1/10, nan,...
    1/11, 1/11, nan,...
    1/12, 1/12, nan,...
    1/13, 1/13];
hold on
plot([0,14], [0.8, 0.8], 'r','linewidth', 2.75)
plot(x,y, 'k', 'linewidth', 2.75)
hold off
xlim([0.5, 13.5])
ylim([0,1.02])
ylabel('Classification Accuracy', 'fontsize', 15)

 set(gca,'FontSize',15)

 legend('Min Exp Accuracy (13 DOF)', 'Chance')
 
xlabel('Degree of Freedom (DOF)', 'fontsize', 15)
set(gca, 'xtick', 1:13)
set(gca, 'xticklabel', 1:13)

set(gcf, 'Units', 'inches')
set(gcf, 'Position', [1 1 6 3])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 3])

%% Electrode Map

%       NaN     2     1     3     4     6     8    10    14   NaN
%        65    66    33    34     7     9    11    12    16    18
%        67    68    35    36     5    17    13    23    20    22
%        69    70    37    38    48    15    19    25    27    24
%        71    72    39    40    42    50    54    21    29    26
%        73    74    41    43    44    46    52    62    31    28
%        75    76    45    47    51    56    58    60    64    30
%        77    78    82    49    53    55    57    59    61    32
%        79    80    84    86    87    89    91    94    63    95
%       NaN    81    83    85    88    90    92    93    96   NaN


elecMap = [NaN     2     1     3     4     6     8    10    14   NaN;...
       65    66    33    34     7     9    11    12    16    18;...
       67    68    35    36     5    17    13    23    20    22;...
       69    70    37    38    48    15    19    25    27    24;...
       71    72    39    40    42    50    54    21    29    26;...
       73    74    41    43    44    46    52    62    31    28;...
       75    76    45    47    51    56    58    60    64    30;...
       77    78    82    49    53    55    57    59    61    32;...
       79    80    84    86    87    89    91    94    63    95;...
      NaN    81    83    85    88    90    92    93    96   NaN];

  
electrodes = [79, 6, 97, 36, 89, 26, 88, 70, 13, 44, 34, 84, 41, 14, 42, 12, 45]';

visMap = ones(size(elecMap))*0.5;

for i = 1:size(visMap,1)
    for j = 1:size(visMap,2)
    
        if isnan(elecMap(i,j))
            visMap(i,j) = 0;
        elseif sum(repmat(elecMap(i,j), length(electrodes), 1) == electrodes)
            visMap(i,j) = 1;
        end % END IF
    end % END FOR
end % END IF

visMap = flipud(visMap);
elecMap2 = flipud(elecMap);

figure(30)
hold on
for i = 1:size(visMap,1)
    for j = 1:size(visMap,2)
    
        patchX = [j-1, j, j, j-1];
        patchY = [i, i, i-1, i-1];
        
        pData = patch(patchX, patchY, 1);
        
        set(pData, 'FaceColor', repmat(visMap(i,j),1,3));
        
        if ~isnan(elecMap(i,j));
            text(j-0.5, i-0.5, sprintf('%d', elecMap2(i,j)), 'HorizontalAlignment', 'Center');
        end % END IF
    end % END FOR
end % END IF
hold off

axis off
axis square