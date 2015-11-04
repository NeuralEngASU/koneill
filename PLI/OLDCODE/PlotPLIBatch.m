<<<<<<< HEAD
pathName = 'E:\data\PLI\DownsampledData';
=======
pathName = 'E:\data\PLI\downsampleddata\';
>>>>>>> origin/master

grids = defgrids;

maxPLIs = zeros(1, length(grids));

<<<<<<< HEAD
% for ii = 1:length(grids)
%     
%     load(sprintf('%s\\g%dmpli_UNR.mat',pathName,ii))
%     
%     tmpPLI = squeeze(mean(p,1));
%     
%     maxPLIs(ii) = max(tmpPLI(:,1));
%     
% end % END FOR grids

maxPLIs = [0.1430    0.1347    0.0962    0.1383    0.1574    0.2154];
maxPLI = ceil(max(maxPLIs)*4)/4;
maxR = 1;
%%
for jj = 1:length(grids)
    
    disp(jj);
    disp('Load Data')
    load(sprintf('%s\\g%dmpli_UNR.mat',pathName,jj))
    disp('Data Loaded')
    PlotPLI2(pathName, grids(jj), p, r, chanpairs, 'space-invader', maxPLI, maxR)
    
    close all
end % END FOR grids

%% Plot Stats Batch


for jj = 1:length(grids)
    disp(jj)
    [sig,p,stats] = check_pli(jj,'E:\data\PLI\DownsampledData');
    save(fullfile(pathName, [grids(jj).subject,'_Stats.mat']), 'sig', 'p', 'stats');
    savefig(fullfile(pathName, sprintf('%s_Stats', grids(jj).subject)));
    print(fullfile(pathName, sprintf('%s_Stats', grids(jj).subject)), '-dpng')
    
    clf;
    disp(sprintf('%d: Done',jj))
end % END FOR

%% Plot space invader signifigance

pathName = 'E:\data\PLI\DownsampledData';

grids = defgrids;

%%
for jj = 1:length(grids)
    
    disp(jj);
    disp('Load Data')
    load(sprintf('%s\\g%dmpli_UNR.mat',pathName,jj))
    load(fullfile(pathName, sprintf('%s_Stats.mat', grids(jj).subject)));
    disp('Data Loaded')
    PlotPLI3(pathName, grids(jj), sig, chanpairs, 'space-invader', mapCol./255)
    
    close all
end % END FOR grids

=======
for ii = 1:length(grids)
    
    load(sprintf('%s\\g%dmpli_UNR.mat',pathName,ii))
    
    maxPLIs(ii) = max(squeeze(mean(p,1)));
    
end % END FOR grids
>>>>>>> origin/master
