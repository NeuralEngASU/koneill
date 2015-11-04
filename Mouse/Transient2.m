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

% EOF