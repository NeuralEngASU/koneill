%Define Grids used in the PLI study

function [gridDef] = GridDef(patNum)

% 2104PP01
gridDef(1).patNum = '2014PP01';
gridDef(1).numChans = 120;
gridDef(1).layoutShortNames = {'BG', 'IHB', 'IHD', 'AI', 'AS', 'PI', 'PS'};
gridDef(1).layoutLongNames = {'Big Grid', 'Interhemispheric Strip B (Double Sided-Side 1)', 'Interhemispheric Strip D (Double Sided-Side 2)',...
                              'Anterior Inferior (Depth)', 'Anterior Superior (Depth)', 'Posterior Inferior (Depth)', 'Posterior Superior (Depth)'}; 
gridDef(1).layout = {reshape(1:64,8,8), reshape(65:80,8,2)', reshape(81:96,8,2)',  [102:-1:97]',...
                     [108:-1:103]', [114:-1:109]',[120:-1:115]'};
gridDef(1).badChan = [];
% gridDef(1).layoutImage = load('2014PP01_LayoutGlobal.png');

% 2014PP02
gridDef(2).patNum = '2014PP02';
gridDef(2).numChans = 94;
gridDef(2).layoutShortNames = {'G', 'AIN', 'PIN', 'OF', 'AT', 'PT'};
gridDef(2).layoutLongNames = {'Grid', 'Anterior Insular Region (Depth)', 'Posterior Insular Region (Depth)',...
                              'Orbitofrontal Region', 'Anterior Temporal Lobe (under)', 'Posterior Temporal Lobe (under)'}; 
gridDef(2).layout = {reshape(1:64,8,8),[70:-1:65]', [76:-1:71]', [77:82], [83:88], [89:94]};
gridDef(2).badChan = [];
% gridDef(2).layoutImage = load('2014PP02_LayoutGlobal.png');

% 2014PP07
gridDef(3).patNum = '2014PP07';
gridDef(3).numChans = 90;
gridDef(3).layoutShortNames = {'LIH', 'RIH', 'AG', 'PG', 'AD', 'LD', 'PD'};
gridDef(3).layoutLongNames = {'Left Interhemispheric Strip', 'Right Interhemispheric Strip', 'Anterior Grid',...
                              'Posterior Grid', 'Anterior Lateral (Depth)', 'Lateral (Depth)', 'Posterior Lateral (Depth)'}; 
gridDef(3).layout = {reshape(1:16,8,2)',reshape(17:32,8,2)', reshape(33:52, 5,4), reshape(53:72, 5,4),...
                     [78:-1:73]', [84:-1:79]',[90:-1:85]'};
gridDef(3).badChan = [];
% gridDef(3).layoutImage = load('2014PP07_LayoutGlobal.png');

end % END FUNCTION

% EOF