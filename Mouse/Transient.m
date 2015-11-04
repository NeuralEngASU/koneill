%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calcium Transient Manipulator
%       Kevin O'Neill
%       2013-12-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Glut v TTX
% Glut v DHPG
% Glut v GABA
% Glut v GABA v TTX


%% Import Data

experiment = 'DHPG';


switch experiment
    case 'Glut'

        % Defining sheet and data range of the column titles
        xlRange = 'A1:BY1';  % Glutamate
        sheet = 1;
        
        % Only grab the column title data
        [~, columnTitles, ~] = xlsread('Glutamate master file.xlsx', sheet, xlRange); % Glutamate
        
        % Reads in the data from the excel spreadsheet
        xlData = xlsread('Glutamate master file.xlsx'); % Glut
        
    case 'TTX'
        % Defining sheet and data range of the column titles
        xlRange = 'A1:Q1';  % Glutamate
        sheet = 1;
        
        [~, columnTitles, ~] = xlsread('Glutamate + TTX master file.xlsx', sheet, xlRange); % TTX
        xlData = xlsread('Glutamate + TTX master file.xlsx'); % TTX
    case 'DHPG'
        % Defining sheet and data range of the column titles
%         xlRange = 'A1:S1';  % DHPG Old
        xlRange = 'A1:W1';  % DHPG new
        sheet = 1;
        
        % Only grab the column title data
        [~, columnTitles, ~] = xlsread('DHPG master file new.xlsx', sheet, xlRange); % DHPG
        
        % Reads in the data from the excel spreadsheet
        xlData = xlsread('DHPG master file new.xlsx'); % DHPG
% return;
    case 'NMDA'
        % Defining sheet and data range of the column titles
        xlRange = 'A1:AA1';  % NMDA
        sheet = 1;
        
        % Only grab the column title data
        [~, columnTitles, ~] = xlsread('NMDA Master.xlsx', sheet, xlRange); % NMDA
        
        % Reads in the data from the excel spreadsheet
        xlData = xlsread('NMDA Master.xlsx'); % NMDA
        
    case 'GABA'
        % Defining sheet and data range of the column titles
        xlRange = 'A1:AL1';  % GAMA
        sheet = 1;
        
        % Only grab the column title data
        [~, columnTitles, ~] = xlsread('GABA Master.xlsx', sheet, xlRange); % NMDA
        
        % Reads in the data from the excel spreadsheet
        xlData = xlsread('GABA Master.xlsx'); % NMDA
        
    case 'GlutNiFed'
        % Defining sheet and data range of the column titles
        xlRange = 'A1:P1';  % GAMA
        sheet = 1;
        
        % Only grab the column title data
        [~, columnTitles, ~] = xlsread('Glut_NiFed Master.xlsx', sheet, xlRange); % NMDA
        
        % Reads in the data from the excel spreadsheet
        xlData = xlsread('Glut_NiFed Master.xlsx'); % NMDA
        
    otherwise
end





% Sets up the string to match in the regular expresion
expr = 'Average';

% Finds where the defined expression occurs in the cell
resultStr = regexp(columnTitles, expr);

% Converts empty cells to 0
for i = 1:length(resultStr)
    if isempty(resultStr{i})
        resultStr{i} = 0;
    end % END IF
end % END FOR

% Converts cell to a logical matrix. This matrix contains the indicies of
% the average control columns
averageIdx = cell2mat(resultStr);
averageIdx = logical(averageIdx(4:end));

time = xlData(:,3); % Removes the time column
xlData = xlData(:, 4:end);% Ignore the first 3 columns as they do not contain intensity measures


%% DeltaF/F

% Find the indicies of the average control columns
idxNum = 1:length(averageIdx);
idxNum = idxNum(averageIdx);

% Initialize a new array of the same size as xlData
deltaF = nan(size(xlData));

% For all but the last average control column, divide the following columns
% by the control to get deltaF/F
for i = 1:length(idxNum)-1
    
    deltaF(:, idxNum(i):idxNum(i+1)-1) = xlData(:,idxNum(i):idxNum(i+1)-1) ./ repmat(xlData(:,idxNum(i)), 1, idxNum(i+1)-idxNum(i));
    
end % END FOR

% Get the deltaF/F values for the last experiment
deltaF(:, idxNum(end):end) = xlData(:,idxNum(end):end) ./ repmat(xlData(:,idxNum(end)), 1, length(averageIdx)-idxNum(end)+1);

%% Flip transients the right direction

% Use a rough slope determine if the transients should be flipped
for i = 1:length(averageIdx)
    
    % Find the index for the last darta point in a column
    tmpIdx = not(isnan(deltaF(:,i)));
    tmpIdx = sum(tmpIdx);
    
    % pseudo slope between teh first point and the last point
    pseudoSlope = deltaF(tmpIdx,i) - deltaF(1,i);
    
    % Iff the slope is negative, flip the transient
    % Example: [4,3,2,1] -> [4,5,6,7]
    if pseudoSlope < 0;     
        deltaF(:,i) = (-1*deltaF(:,i)) + 2*deltaF(1,i);
    end % END IF
end % END IF


%% Find Start Times
% This could not be done automatically, so I had to do it manually.
% There is an automatic method, but it would not work for all cases.

% Super smooth the data set.
% derp = TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(deltaF))));
derp = (TransientSmooth(deltaF));

% % This is an expeirmental section
% slope = TransientSlope(derp);
% 
% accel = TransientSlope(slope);
% 
% % hold on
% % % plot(deltaF(:,2))
% % plot(slope(:,:)*100)
% % % plot(accel(:,1:3)*200-1, 'r')
% plot(derp)
% % hold off
% 
% [pks, locs] = findpeaks(slope(:,2)*100, 'MINPEAKHEIGHT', 0.1);

% Plot all of the transients and manually determine the start time
% for i = 1:length(averageIdx)
%     
%     figure(i)
% %     plot(deltaF(:,i))
%     plot(derp(:,i))
%     title(sprintf('Idx: %d', i));
%     
% end % END FOR

switch experiment
    case 'Glut'
        
        timeIdx   = [74, 73, 69, 67, 66, 67, 66, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 39, 38, 37, 36, 31, 30, 27, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 13, 12, 9,  8,  7,  6,  5,  4 , 3,  2, 1];
        timeStart = [20,  0,  0, 64,  0,  4 , 0,  0, 36,  0,  5,  5,  5,  5,  5,  0, 15, 15, 15, 15, 15,  0, 21, 21, 21, 21, 21,  0, 55,  0, 35, 35,  0, 85,  0,  0,  0,  0,  0,  0,  0,  0, 36, 36, 36, 36,  0, 52,  0, 0,  0, 0, 36, 36, 36, 36, 36, 0];
        
        toFlip = [74];
        exclude = [72, 71, 70, 68, 65, 64, 63, 40, 35, 34, 33, 32, 29, 28, 26, 14, 11, 10];
        
    case 'TTX'
        
        timeIdx   = [14,13,12,11,10];
        timeStart = [42,42,42,42, 0];
        
        toFlip = [];
        exclude = [9, 8, 7, 6, 5, 4, 3, 2, 1];
        
    case 'DHPG'
                
        timeIdx   = [ ];
        timeStart = [ ];
        
        toFlip = [];
        exclude = [16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1];
        
    case 'NMDA'
    otherwise
end

% return;

%% Re-flip selected microglia
% for i = toFlip
%     
%     deltaF(:,toFlip) = (-1*deltaF(:,i)) + 2*deltaF(1,i);
% 
% end % END FOR
%% Normalize to 0

% Create a new matrix the same size as deltaF
normData = nan(size(deltaF));

% For every index as found in the 'Find Start Time' section, normalize to
% zero
for i = 1:length(timeIdx)
    
    % Column index
    c = timeIdx(i);
    
    if timeStart(i) == 0 % indicies cannot start at 0 in MATLAB
        r = 1; % Row index
        timeStart(i) = 1;
    elseif timeStart(i) > 0
        r = 1:timeStart(i)-1;
    else
        disp('Something went horribly wrong.')
    end % END IF
     
    % Normalization = data/baseline - 1
    normData(:,c) = deltaF(:,c) ./ nanmean(deltaF(r,c)) - 1;
    
end % END FOR

% Normalization for excluded ROIs
for i = 1:length(exclude)
    
    % Column index
    c = exclude(i);
    r = 1;
     
    % Normalization = data/baseline - 1
    normData(:,c) = deltaF(:,c) ./ nanmean(deltaF(r,c)) - 1;
    
end % END FOR

%% Adjust start times

% Find the minimim and maximum start times
maxTime = max(timeStart);
minTime = min(timeStart);

% Make dure the minimum start time is at least 0
if minTime == 0
    minTime = 1;
end % END IF

% Calculate the difference between the min and max start time
diffTime = maxTime - minTime;

% Make a new matrix the same size as normData, but with an additional
% setion of NaNs at the bottom
normData2 = [nan(size(normData)); nan(diffTime, size(normData, 2))];

% Adjusts the time for the non-excluded ROIs
for i = 1:length(timeIdx)
    
    c = timeIdx(i);
    
    if sum(exclude == c) % Set excluded ROIs to NaN
        normData(:,c) = nan(size(normData,1), 1);
    elseif not(nansum(normData(:,c))) % I forgot what this did but it's important
        normData(:,c) = nan(size(normData,1), 1);
    else
        %                [          NaN Packing,             Data,           NaN Packing             ]
        normData2(:,c) = [nan(maxTime - timeStart(i),1); normData(:,c); nan(timeStart(i) - minTime,1)];
        
    end % END IF
    
end % END FOR

%% Average Analysis

% testData = nanmean(normData2, 2);

% Mean across all of the data (average controls and excluded ROIs are
% already removed)
testData = nanmean(normData2(:,:), 2);

% Find the standard error across the rows (stdev/sqrt(n))
testErr  = nanstd(normData2(:,:), 1, 2);
testErr  = testErr ./ sqrt(sum(not(isnan(normData2(100, :)))));

x = 1:length(normData2(:,1));
xx = [x, fliplr(x)];

x  =  x * 0.99 - maxTime * 0.99; % Time with first uncaging at 0
xx = xx * 0.99 - maxTime * 0.99; % Time with first uncaging at 0


stestData = TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(testData))))));
snormData = TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(normData2))))));

stestErr  =  nanstd(snormData(:,:), 1, 2);
stestErr = stestErr ./ sqrt(sum(not(isnan(snormData(100, :)))));

patchData = [stestData' + stestErr', fliplr(stestData' - stestErr')];

% x = [nan(1, maxTime - 36), time', nan(1, 36 - minTime)];
% xx = [x, fliplr(x)];
xx = xx(not(isnan(patchData)));

patchData = patchData(not(isnan(patchData)));

% tmpTime = [0:0.99:(floor(sum(not(isnan(patchData)))/2)-1)*.99];
% xx = [tmpTime, fliplr(tmpTime)];

% tmpTime = [nan(1, 35), tmpTime];
% tmpTime = [tmpTime, nan(1, length(stestData)-length(tmpTime))];

figure(10)

% plot(testData)
hold on
pData = patch(xx/60, patchData,1);
plot(x/60, stestData, 'r', 'linewidth', 2.75)
hold off

set(pData, 'FaceAlpha', 0.25)
set(pData, 'FaceColor', 'r')
set(pData, 'EdgeColor', 'none')

title(sprintf('Microglia (24May) Averaged Data\nGlutemate Uncaging'))
ylabel('%\DeltaF/F')
xlabel('Time, min')

ylim([-1, 5])
xlim([floor(x(1)/60), ceil(x(end)/60)])

set(gca, 'YTick', [-1, 0, 1, 2, 3, 4, 5, 6])
set(gca, 'YTickLabel', [-100, 0, 100, 200, 300, 400, 500, 600]);

set(gcf, 'Units', 'inches')
set(gcf, 'Position', [1 1 5 2.5])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 2.5])



print('-dpng', ['C:\CodeRepo\Lab\Mouse\Mouse\Figures\', experiment ,'_AvgStdErr_Patch_24May.png'], '-r100');

close gcf

%% Error Bar Plot All Data

% Data manipulation
testData = nanmean(normData2(:,:), 2);

stestData = TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(testData))))));
snormData = TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(normData2))))));

stestErr  =  nanstd(snormData(:,:), 1, 2);
stestErr = stestErr ./ sqrt(sum(not(isnan(snormData(100, :)))));

% Plot
figure(20)
clear('gcf')
hold on
plot(x(1:25:end)/60, stestData(1:25:end), '-or', 'MarkerSize',3, 'MarkerFaceColor', 'r');

for i = 1:25:length(x)
    
    plot([x(i), x(i)]/60, [stestData(i) + stestErr(i), stestData(i) - stestErr(i)], '-r') 
    
end

hold off

title(sprintf('Averaged Glutemate Uncaging Data\nAll Data: n = %d', sum(not(isnan(normData2(100, :)))) ))
ylabel('%\DeltaF/F')
xlabel('Time, min')

ylim([-1, 5])
xlim([floor(x(1)/60), ceil(x(i)/60)])

set(gca, 'YTick', [-1, 0, 1, 2, 3, 4, 5, 6])
set(gca, 'YTickLabel', [-100, 0, 100, 200, 300, 400, 500, 600]);

set(gcf, 'Units', 'inches')
set(gcf, 'Position', [1 1 5 2.5])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 2.5])


print('-dpng', ['C:\CodeRepo\Lab\Mouse\Mouse\Figures\', experiment ,'_AvgStdErr_ErrorBar_All.png'], '-r100');

close gcf
% errorbar(x(1:20:end)/60, stestData(1:20:end), stestErr(1:20:end), 'o-r', 'MarkerFaceColor', 'r', 'MarkerSize',3)

%% Errorbar Plot: 24May 

% Data manipulation


switch experiment
    case 'Glut'
        range = 2:6;
    case 'TTX'
        range = 11:14;
    otherwise
end

testData = nanmean(normData2(:,range), 2);

stestData = TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(testData))))));
snormData = TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(normData2))))));

stestErr  =  nanstd(snormData(:,range), 1, 2);
stestErr = stestErr ./ sqrt(sum(not(isnan(snormData(100, range)))));

% Find Transients
slope = TransientSlope(snormData(maxTime:end,:));

tmpX = (1:length(normData(:,2))) - 36;
% tmpX = repmat(tmpX', 1, 4);
% plot(tmpX, normData(:,[2,3,4,5,6]))

% plot(x/60, snormData(:,[2,3,4,5,6]))

[pks, locs] = findpeaks(slope(:,2)*100, 'MINPEAKHEIGHT', 2.0);

distStim = locs(1) - maxTime;

% Plot
figure(21)
clear('gcf')
hold on
plot(x(1:25:end)/60, stestData(1:25:end), '-or', 'MarkerSize',3, 'MarkerFaceColor', 'r');

% Plot errorbars
for i = 1:25:length(x)
    plot([x(i), x(i)]/60, [stestData(i) + stestErr(i), stestData(i) - stestErr(i)], '-r') 
end % END FOR

for j = 1:length(locs)
    
    j = locs(j);
        
    plot([x(j)-distStim, x(j)+distStim]./60, [6, 6], 'k', 'LineWidth', 5)
    
end % END FOR

hold off

title(sprintf('Averaged Glutemate Uncaging Data\n24May Data: n = %d', sum(not(isnan(normData2(100, range)))) ))
ylabel('%\DeltaF/F')
xlabel('Time, min')

ylim([-1, 6])
xlim([floor(x(1)/60), ceil(x(i)/60)])

set(gca, 'YTick', [-1, 0, 1, 2, 3, 4, 5, 6])
set(gca, 'YTickLabel', [-100, 0, 100, 200, 300, 400, 500, 600]);
text(0, 5.5, 'Uncaging')

set(gcf, 'Units', 'inches')
set(gcf, 'Position', [1 1 5 2.5])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 2.5])


print('-dpng', ['C:\CodeRepo\Lab\Mouse\Mouse\Figures\AvgStdErr_ErrorBar_24May.png'], '-r100');



%% Baseline Jump Analysis

trialTime = maxTime;

baseline = nanmean(normData(1:36,range), 1);

transMat = nan(length(locs), length(range));
errMat = nan(length(locs), length(range));

distStim = 3;


for transient = 1:length(locs);
    
    idx = locs(transient) + distStim;
    transMat(transient,:) = snormData(idx, range);   
    
end % END FOR

errMat = nanstd(transMat, 1, 2);
errMat = errMat ./ sqrt(sum(not(isnan(snormData(100, range)))));

% For use if we want to re-normalize
% baseMat = repmat(baseline, length(locs), 1);

transMean = nanmean(transMat,2);

figure(30)
hold on
plot(1:length(locs), transMean, 'or', 'MarkerSize', 3, 'MarkerFaceColor', 'r')
for i = 1:length(locs)
    plot([i, i], [transMean(i) + errMat(i), transMean(i) - errMat(i)], '-r') 
end % END FOR
hold off

title(sprintf('Uncaged Intensity Peak\nn = %d', sum(not(isnan(normData2(100, range)))) ))
ylabel('%\DeltaF/F')
xlabel('Uncage Event')

ylim([0, 6])
xlim([1, length(transMean)])

set(gca, 'YTick', [-1, 0, 1, 2, 3, 4, 5, 6])
set(gca, 'YTickLabel', [-100, 0, 100, 200, 300, 400, 500, 600]);

set(gca, 'XTick', [1:8])
set(gca, 'XTickLabel', [1:8]);

set(gcf, 'Units', 'inches')
set(gcf, 'Position', [1 1 3 2.5])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3 2.5])

print('-dpng', ['C:\CodeRepo\Lab\Mouse\Mouse\Figures\', experiment, '_AvgStdErr_ErrorBar_UncageStep_24May.png'], '-r100');

%%


% plot(slope(:, timeIdx(2:end)))
% plot(slope(:, 8))
% 
% trial = 2;
% hold on
% plot(normData2(:,trial))
% plot(slope(:, trial)*20, 'r')
% hold off

% xlim([0, 500])
% 
% derp = TransientSlope(TransientSmooth(TransientSmooth(TransientSmooth(normData2))));

% plot(derp(:, trial))

% [pks, locs] = findpeaks(slope(:, trial), 'MINPEAKHEIGHT', max(slope(:,8))*.25)

% accel = TransientSlope(slope);

% plot(accel(:,2))

% ylim([-.02, .02])

% plot(snormData(:, timeIdx(1:end)))
% 

slope = TransientSlope(snormData);

switch experiment
    case 'Glut'
        peakMat   = nan(8, length(averageIdx));
        transMat2 = nan(8, length(averageIdx));
        maxTrans = 8;
    case 'TTX'
        peakMat   = nan(8, length(averageIdx));
        transMat2 = nan(8, length(averageIdx));
        maxTrans = 8;
    otherwise
end

sortedPks = [];
sortedLocs = [];
pks = [];
locs = [];

for i = 1:length(timeIdx)
    
    trial = timeIdx(i);
    
    if ~averageIdx(trial)
        
        timeOffset = maxTime-10;
        
        [pks, locs] = findpeaks(slope(timeOffset:end, trial));
        
        locs = locs + timeOffset;
        
        [sortedPks sortedLocs] = sort(pks, 'descend');
        
        if length(sortedPks) < maxTrans
            disp(sprintf('Issue with trial: %d', trial));
        else
            sortedPks = sortedPks(1:maxTrans);
            sortedLocs = locs(sortedLocs(1:maxTrans));
%             sortedLocs = sortedLocs(1:maxTrans);
            
            origLocs = sort(sortedLocs, 'ascend');
            
            peakMat(:, trial) = origLocs(:);
            
        end % END IF
    end % END IF
end % END IF

disp('pass')

peakMat = peakMat + 3;

for i = 1:maxTrans
    for j = 1:length(timeIdx)
        
        j = timeIdx(j);
        
        if ~averageIdx(j)
            transMat2(i,j) = snormData(peakMat(i,j), j);
        end % END IF
        
    end % END IF
end % END IF


errMat2 = nanstd(transMat2, 1, 2);
errMat2 = errMat2 ./ sqrt(sum(not(isnan(transMat2(1,:)))));

% For use if we want to re-normalize
% baseMat = repmat(baseline, length(locs), 1);

transMean2 = nanmean(transMat2,2);


switch experiment
    case 'Glut'
        glutTransMean = transMean2;
        glutErrMat = errMat2;
        glutNumTrials = sum(not(isnan(transMat2(1, :))));
        
        save('GlutData', 'glutTransMean')
        save('GlutData', 'glutErrMat', '-append')
        save('GlutData', 'glutNumTrials', '-append')
        
    case 'TTX'
        ttxTransMean = transMean2;
        ttxErrMat = errMat2;
        ttxNumTrials = sum(not(isnan(transMat2(1, :))));
    otherwise
end

figure(32)
hold on
plot(1:8, transMean2, 'or', 'MarkerSize', 3, 'MarkerFaceColor', 'r')
for i = 1:8
    plot([i, i], [transMean2(i) + errMat2(i), transMean2(i) - errMat2(i)], '-r') 
end % END FOR
hold off

title(sprintf('Uncaged Intensity Peak\nn = %d', sum(not(isnan(transMat2(1, :)))) ))
ylabel('%\DeltaF/F')
xlabel('Uncage Event')

ylim([0, 5])
xlim([1, length(transMean2)])

set(gca, 'YTick', [-1, 0, 1, 2, 3, 4, 5, 6])
set(gca, 'YTickLabel', [-100, 0, 100, 200, 300, 400, 500, 600]);

set(gca, 'XTick', [1:8])
set(gca, 'XTickLabel', [1:8]);

set(gcf, 'Units', 'inches')
set(gcf, 'Position', [1 1 3 2.5])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3 2.5])

print('-dpng', ['C:\CodeRepo\Lab\Mouse\Mouse\Figures\', experiment ,'_AvgStdErr_ErrorBar_UncageStep_All.png'], '-r100');



%%
% % return;

if exist('glutTransMean') && exist('ttxTransMean')
    
    figure(33)
    hold on
    
    % plots for legend

    plot(1:7, glutTransMean(1:end-1), '-or', 'MarkerSize', 3, 'MarkerFaceColor', 'r')
    plot(1:7, ttxTransMean, '-og', 'MarkerSize', 3, 'MarkerFaceColor', 'g')

    
    
        % Plot Glut
    plot(1:7, glutTransMean(1:end-1), '-or', 'MarkerSize', 3, 'MarkerFaceColor', 'r')
    for i = 1:7
        plot([i, i], [glutTransMean(i) + glutErrMat(i), glutTransMean(i) - glutErrMat(i)], '-r')
    end % END FOR
    
    % Plot Glut + TTX
    plot(1:7, ttxTransMean, '-og', 'MarkerSize', 3, 'MarkerFaceColor', 'g')
    for i = 1:7
        plot([i, i], [ttxTransMean(i) + ttxErrMat(i), ttxTransMean(i) - ttxErrMat(i)], '-g')
    end % END FOR
    

    

    
    hold off
    
    title(sprintf('Uncaged Intensity Peak Comparison\n(Glut: n = %d) (Glut+TTX: n = %d)', glutNumTrials, ttxNumTrials ))
    ylabel('%\DeltaF/F')
    xlabel('Uncage Event')
    
    legend({'Glut', 'Glut + TTX'}, 'Location', 'Northwest')
    
    ylim([0, 5])
    xlim([1, 7])
    
    set(gca, 'YTick', [-1, 0, 1, 2, 3, 4, 5, 6])
    set(gca, 'YTickLabel', [-100, 0, 100, 200, 300, 400, 500, 600]);
    
    set(gca, 'XTick', [1:8])
    set(gca, 'XTickLabel', [1:8]);
    
    set(gcf, 'Units', 'inches')
    set(gcf, 'Position', [1 1 3 2.5])
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3 2.5])
    
    set(gca,'TickDir','out')
    
    print('-dpng', ['C:\CodeRepo\Lab\Mouse\Mouse\Figures\GlutvsTTX_AvgStdErr_ErrorBar_UncageStep_All.png'], '-r100');
    
end



% [h,p] = ttest(glutTransMean(1:end-1), ttxTransMean);

%%

% return;

%%

if isequal(experiment, 'TTX')

ttxIdx = [ NaN 188 188 188 188 NaN 111 111,111 NaN  60  60  60  60;...
           NaN 222 222 222 222 NaN 154 154 154 NaN  72  72  72  72;...
           NaN 273 273 273 273 NaN 211 211 211 NaN  97  97  97  97;...
           NaN 339 339 339 339 NaN 283 283 283 NaN 127 127 127 127;...
           NaN 389 389 389 389 NaN 326 326 326 NaN 140 140 140 140;...
           NaN 448 448 448 448 NaN 404 404 404 NaN 161 161 161 161;...
           NaN 507 507 507 507 NaN 441 441 441 NaN 186 186 186 186];
       
ttxStart = [1 140 140 140 140 1 87 87 87 1 42 42 42 42];

ttxIdx(:,7) = nan(7,1);


normData3 = nan(size(deltaF));

for i = 1:size(deltaF, 2)
    
    baseline = nanmean(deltaF(1:ttxStart(i), i));
    
    normData3(:,i) = deltaF(:,i) ./ baseline - 1;
    
end % END FOR

snormData = TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(normData3))))));

jumpData = ttxIdx;

for i = 1:7
    for j = 1:size(deltaF, 2)
        
        if ~isnan(ttxIdx(i,j))
            jumpData(i,j) = snormData(ttxIdx(i,j),j);
        end % END IF
        
    end % END FOR
end % END FOR

ttxTransMean = nanmean(jumpData(:,:), 2);

ttxErrMat  =  nanstd(jumpData(:,:), 1, 2);
ttxErrMat = ttxErrMat ./ sqrt(sum(not(isnan(jumpData(1,:)))));

ttxNumTrials = sum(not(isnan(jumpData(1,:))));

save('TTXData', 'ttxTransMean')
save('TTXData', 'ttxErrMat', '-append')
save('TTXData', 'ttxNumTrials', '-append')
save('TTXData', 'ttxIdx', '-append')
save('TTXData', 'ttxStart', '-append')


end


%%

if isequal(experiment, 'DHPG')

dhpgIdx = [NaN,60,60,60,NaN,29,NaN,NaN,NaN,NaN,NaN,12,12,NaN,18,18,NaN,18,...
    18,18;NaN,83,83,83,NaN,44,NaN,NaN,NaN,NaN,NaN,37,37,NaN,49,49,NaN,34,...
    34,34;NaN,105,105,105,NaN,58,NaN,NaN,NaN,NaN,NaN,51,51,NaN,70,70,NaN,...
    61,61,61;NaN,129,129,129,NaN,74,NaN,NaN,NaN,NaN,NaN,65,65,NaN,94,94,...
    NaN,83,83,83;NaN,153,153,153,NaN,88,NaN,NaN,NaN,NaN,NaN,88,88,NaN,116,...
    116,NaN,111,111,111;NaN,178,178,178,NaN,104,NaN,NaN,NaN,NaN,NaN,110,...
    110,NaN,129,129,NaN,138,138,138;NaN,202,202,202,NaN,122,NaN,NaN,NaN,...
    NaN,NaN,137,137,NaN,147,147,NaN,156,156,156];
       
dhpgStart = [1	12	12	12	1	1	1	7	7	7	1	1	1	1	1	1	1	6	6	6];

normData3 = nan(size(deltaF));
for i = 1:size(deltaF, 2)
    
    baseline = nanmean(deltaF(1:dhpgStart(i), i));
    
    normData3(:,i) = deltaF(:,i) ./ baseline - 1;
    
end % END FOR

snormData = TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(normData3))))));

jumpData = dhpgIdx;

for i = 1:7
    for j = 1:size(deltaF, 2)
        
        if ~isnan(dhpgIdx(i,j))
            jumpData(i,j) = snormData(dhpgIdx(i,j),j);
        end % END IF
        
    end % END FOR
end % END FOR

dhpgTransMean = nanmean(jumpData(:,:), 2);

dhpgErrMat  =  nanstd(jumpData(:,:), 1, 2);
dhpgErrMat = dhpgErrMat ./ sqrt(sum(not(isnan(jumpData(1,:)))));

dhpgNumTrials = sum(not(isnan(jumpData(1,:))));

save('DHPGData', 'dhpgTransMean')
save('DHPGData', 'dhpgErrMat', '-append')
save('DHPGData', 'dhpgNumTrials', '-append')
save('DHPGData', 'dhpgIdx', '-append')
save('DHPGData', 'dhpgStart', '-append')


end

%% GABA

if isequal(experiment, 'GABA')

gabaIdx = [NaN,14,14,14,14,NaN,10,10,10,10,10,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,...
    13,13,13,13,NaN,16,16,16,16,NaN,27,NaN,6,6,6;NaN,37,37,37,37,NaN,19,19,19,19,...
    19,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,24,24,24,24,NaN,44,44,44,44,NaN,40,...
    NaN,29,29,29;NaN,53,53,53,53,NaN,41,41,41,41,41,NaN,NaN,NaN,NaN,NaN,NaN,NaN,...
    NaN,NaN,41,41,41,41,NaN,62,62,62,62,NaN,50,NaN,56,56,58;NaN,73,73,73,73,NaN,...
    53,53,53,53,53,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,72,72,72,72,NaN,78,78,78,...
    78,NaN,73,NaN,68,68,66;NaN,95,95,95,95,NaN,84,84,84,84,84,NaN,NaN,NaN,NaN,...
    NaN,NaN,NaN,NaN,NaN,91,91,91,91,NaN,110,110,110,110,NaN,88,NaN,79,79,79;...
    NaN,122,122,122,122,NaN,99,99,99,99,99,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,...
    NaN,106,106,106,106,NaN,125,125,125,125,NaN,117,NaN,100,100,100;NaN,147,...
    147,147,147,NaN,113,113,113,113,113,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,...
    121,121,121,121,NaN,157,157,157,157,NaN,139,NaN,111,111,111];
       
gabaStart = [1	8	8	8	8	1	8	8	8	8	8	1	1	1	1	1	1	1	1	1	7	7	7	7	1	6	6	6	6	1	18	1	3	3	3];

gabaIdx(:,7) = nan(7,1);

normData3 = nan(size(deltaF));
for i = 1:size(deltaF, 2)
    
    baseline = nanmean(deltaF(1:gabaStart(i), i));
    
    normData3(:,i) = deltaF(:,i) ./ baseline - 1;
    
end % END FOR

snormData = TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(TransientSmooth(normData3))))));

jumpData = gabaIdx;

for i = 1:7
    for j = 1:size(deltaF, 2)
        
        if ~isnan(gabaIdx(i,j))
            jumpData(i,j) = snormData(gabaIdx(i,j),j);
        end % END IF
        
    end % END FOR
end % END FOR

gabaTransMean = nanmean(jumpData(:,:), 2);

gabaErrMat  =  nanstd(jumpData(:,:), 1, 2);
gabaErrMat = gabaErrMat ./ sqrt(sum(not(isnan(jumpData(1,:)))));

gabaNumTrials = sum(not(isnan(jumpData(1,:))));

save('GABAData', 'gabaTransMean')
save('GABAData', 'gabaErrMat', '-append')
save('GABAData', 'gabaNumTrials', '-append')
save('GABAData', 'gabaIdx', '-append')
save('GABAData', 'gabaStart', '-append')


end


%%
% for i = 1:25
%     figure(i)
%     
%     plot(x, snormData(:,i+1))
%     title(sprintf('Microglia (24May) ROI %d Data\nGlutemate Uncaging', i))
%     
%     ylim([-0.5, 6])
%     xlim([0, 922])
% 
%     set(gca, 'YTick', [0, 1, 2, 3, 4, 5, 6])
%     set(gca, 'YTickLabel', [0, 100, 200, 300, 400, 500, 600]);
%     
%     ylabel('%\DeltaF/F')
%     xlabel(' Time, sec')
% 
%     print('-dpng', ['C:\CodeRepo\Lab\Mouse\Mouse\Figures\ROI', num2str(i),'.png'], '-r100');
% 
%     
% end % END IF

% EOF