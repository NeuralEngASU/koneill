%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	MouseCorrCoef
%		Kevin O'Neill
%		Mario Capecchi Lab
%		20140312
%		v0.1
%		PI: Naveen Nagarajan
%
%	Inputs:
%       CSCData: Data loaded from the expYYYY-MM-DD_hh-mm-ss_CSC.mat file
%											(ns2mat output)
%
%       saveDate:  	Saves plots or regenerated plots to the hard drive. Data will be saved as
%					a .fig, .png, and a .pdf. File names will follow the convention:
%						CorrCoef_YYYY-MM-DD_hh-mm-ss.ext
%					saveData = 1: Save data
%					saveData = 0: Do not save data
%
%       stimData:   Data containing the laser stimulaton if it is needed.
%
%		timeBounds:	Plot between the times specified (inclusive) around an event. Must be positive if events is empty.
%					Example:   [100, 1000] % seconds, no events
%					Example:   [-25, 25]   % Seconds, plot from -25 before and 25 after event
%
%		events:		If events is not empty, plot around event (such as stim turned on) and use timeBounds
%					as the plot boundaries. Events are a vector of timeStamps referencing when stimulation
%					turned on/off (give only one side of the simulation. timeBounds will tell the plotter
%					what time section around the event to plot. Units are in seconds
%					Example: trueEvents = [25, 50, 75, 100, 125, 150, ...]   on/off, on/off, on/off
%						     events = [ 25,  75, 125]   stim on events only
%							 events = [ 50, 100, 150]   stim off events only
%							 timeBounds = [-50, 100]    plot time from [event-50] to [event+100]
%
%       nameMod:    Name modifier. Appends the given string to the end of
%                   the saved file name. Can be a cell if multiple names
%                   will be used. Length of cell must be equal to the
%                   number of events in this case.
%
%       calcOption: The default is option 1 where the code will preform a
%                   corrcoef for each separate event. Option 2 is the
%                   combined corrcoef for all events. In either case, the
%                   complete corrcoef matrix will be output.
%
%       plotOption: Changes the corrcoef plotting style. A value of 1 will
%                   preform the default 'array within an array' grid. An
%                   option of 2 will plot Rij where ij are the two channels
%                   being compared.
%
%	Outputs:
%       1) Figures may be saved by using the saveData flag.
%
%       2) R: Numeric outputs are given in planes formated as Rij, where ij
%             are the two channels being compared. Each plane will be a 
%             separate event.
%
%	To Use:
%		Run function. All inputs must be given.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ R ] = MouseCorrCoef( CSCData, saveData, stimData, timeBounds, events, nameMod, calcOption, plotOption)


% Ensures that there is at least one arbitrary event
if ~isempty(events)
    numEvents = size(events);
else
    numEvents = 1;
    events = 0;
end % END IF

% Sets the number of channels to be used.
if stimData == 1
    numChans = 17;
elseif stimData == 0;
    numChans = 16;
end % END IF

% Preallocates the corrcoef matrix
R = zeros(numChans, numChans, numEvents);

% Loops over the number of events given. If no event was given, then loop
% over the single time given.
for i = 1:numEvents
    startTime = events(1) + timeBounds(1);	% Seconds
    endTime   = events(1) + timeBounds(2);	% Seconds
    
    startTime = startTime * 10e5;			% micro seconds
    endTime   = endTime   * 10e5;           % micro seconds
    
    timeStamps = CSCData.timeStamps; % vector containing the timestamps for all channels
    timeIdx = timeStamps > startTime & timeStamps < endTime;
    
    % Collects the channel data into a column matrix
    for chan = 1:numChans
        
        tempChanName = ['CSC',num2str(chan)];
        xData(:,i) = h5read(CSCData.hdf5FileName, ['/', tempChanName], timeIdx(1), timeIdx(end) - timeIdx(1));
        
    end % END FOR
    
    R(:,:,i) = corrcoef(xData);
    
end % END FOR

%% Plot CorrCoef

if calcOption == 2
    R2 = average(R, 3);
    
    Rtemp = R;
    R = R2;
end % END IF

%% Plot Option 1

if plotOption == 1;
    
    figure(1)
    clf
    
    colorDensity = 128; % Sets the number of colors desired.
    mapColor = jet(colorDensity); % Contains the RGB valuse for the jet color map.
    
    sqPchCoord = GenPatchCoord1(numChans);
    

    
    for i = 1:numChans
        for j = 1:numChans
            
            smallPatch = patch(sqPchCoord(:,j,1), sqPchCoord(:,i,2), 1);
            
            set(smallPatch, 'FaceColor', mapColor(ceil((R(i,j)+1)/2 * colorDensity),:))
            set(smallPatch, 'EdgeColor', 'none')
            
        end % END FOR
    end % END FOR
    
 % Plot Option 2
    
elseif plotOption == 2
    % Electrode 'physical' arrangement:
    %  1  2  3  4
    %  5  6  7  8
    %  9 10 11 12
    % 13 14 15 16
    %
    % Do we want a 10x10 UEA grid image with non-active electrodes in black?
    figure(1)
    clf
    
    if numChans == 17
        numSquare = 16*16; % Each electrode will have 16 plots associated with it.
    else
        numSquare = numChans*numChans;
    end % END IF
    
    % colormap(jet); % Sets the current color map to jet
    colorDensity = 128; % Sets the number of colors desired.
    mapColor = jet(colorDensity); % Contains the RGB valuse for the jet color map.
    
    sqPchCoord = GenPatchCoord2(numChans);
    
    hold on
    for i = 1:16
        
        for j=1:16
            
            squarePatchxx = [sqPchCoord(i,1,1), sqPchCoord(i,1,1)+1, sqPchCoord(i,1,1)+1, sqPchCoord(i,1,1)];
            squarePatchyy = [sqPchCoord(i,1,2), sqPchCoord(i,1,2), sqPchCoord(i,1,2)-1, sqPchCoord(i,1,2)-1];
            
            switch j
                case 1
                    squarePatchxx = squarePatchxx;
                    squarePatchyy = squarePatchyy;
                case 2
                    squarePatchxx = squarePatchxx + 1;
                    squarePatchyy = squarePatchyy;
                case 3
                    squarePatchxx = squarePatchxx + 2;
                    squarePatchyy = squarePatchyy;
                case 4
                    squarePatchxx = squarePatchxx + 3;
                    squarePatchyy = squarePatchyy;
                case 5
                    squarePatchxx = squarePatchxx;
                    squarePatchyy = squarePatchyy - 1;
                case 6
                    squarePatchxx = squarePatchxx + 1;
                    squarePatchyy = squarePatchyy - 1;
                case 7
                    squarePatchxx = squarePatchxx + 2;
                    squarePatchyy = squarePatchyy - 1;
                case 8
                    squarePatchxx = squarePatchxx + 3;
                    squarePatchyy = squarePatchyy - 1;
                case 9
                    squarePatchxx = squarePatchxx;
                    squarePatchyy = squarePatchyy - 2;
                case 10
                    squarePatchxx = squarePatchxx + 1;
                    squarePatchyy = squarePatchyy - 2;
                case 11
                    squarePatchxx = squarePatchxx + 2;
                    squarePatchyy = squarePatchyy - 2;
                case 12
                    squarePatchxx = squarePatchxx + 3;
                    squarePatchyy = squarePatchyy - 2;
                case 13
                    squarePatchxx = squarePatchxx;
                    squarePatchyy = squarePatchyy - 3;
                case 14
                    squarePatchxx = squarePatchxx + 1;
                    squarePatchyy = squarePatchyy - 3;
                case 15
                    squarePatchxx = squarePatchxx + 2;
                    squarePatchyy = squarePatchyy - 3;
                case 16
                    squarePatchxx = squarePatchxx + 3;
                    squarePatchyy = squarePatchyy - 3;
                otherwise
            end % END SWITCH
            
            smallPatch = patch(squarePatchxx, squarePatchyy, 1);
            
            set(smallPatch, 'FaceColor', mapColor(round((R(i,j)+1)/2 * colorDensity),:))
            set(smallPatch, 'EdgeColor', 'none')
            
        end % END FOR
    end % END FOR
    
    for i=1:16
        %     squarePatchxx = [sqPchCoord(i,2), sqPchCoord(i,2)+2, sqPchCoord(i,2)+2, sqPchCoord(i,2)];
        %     squarePatchyy = [sqPchCoord(i,1), sqPchCoord(i,1), sqPchCoord(i,1)-2, sqPchCoord(i,1)-2];
        
        squarePatchxx = [sqPchCoord(i,1,1), sqPchCoord(i,1,1)+4, sqPchCoord(i,1,1)+4, sqPchCoord(i,1,1)];
        squarePatchyy = [sqPchCoord(i,1,2), sqPchCoord(i,1,2), sqPchCoord(i,1,2)-4, sqPchCoord(i,1,2)-4];
        
        largePatch = patch(squarePatchxx, squarePatchyy, 1);
        
        set(largePatch, 'FaceColor', 'none')
        set(largePatch, 'EdgeColor', 'k')
    end % END FOR
    
    
else
    fprintf('The plot option is not recognized. Please use one of the recognized options.\n');
    return;
end % END IF

colorbar('EastOutside')
colorbar('YTickLabel',linspace(-1,1,11))


title('Electrode Cross-Correlation Test')
hold off
xlim([0, numChans])
ylim([0, numChans])

set(gca, 'XTick', [0])
set(gca, 'XTickLabel',{''})

set(gca, 'YTick', [0])
set(gca, 'YTickLabel',{''})

set(gca, 'Visible', 'Off')

%%

% set(gcf, 'Units', 'inches')
% set(gcf, 'Position',[0 0 5 4])
% set(gcf, 'PaperUnits','inches','PaperPosition',[0 0 5 4])

shortName = ['CorrCoef_', CSCData.expDate, '_', nameMod{i}];

if saveData
    print('-dpng', ['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.png'], '-r100');
    saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.fig']);
    saveas(gcf,['G:\CodeRepo\MATLAB\Mouse\Analysis\Figures\', shortName, '.pdf']);
end % END IF

R = Rtemp;

end % END FUNCTION

function [ patchCoord ] = GenPatchCoord1( numChans )

patchCoord = zeros(4, numChans, 2);

xVals = zeros(4, numChans, 1);
yVals = zeros(4, numChans, 1);

plotHeight = [numChans:-1:0];
plotWidth  = [0:1:numChans];
for i = 1:numChans
    
    xVals(1,i) = plotWidth(i);
    xVals(2,i) = plotWidth(i+1);
    xVals(3,i) = plotWidth(i+1);
    xVals(4,i) = plotWidth(i);
    
    yVals(1,i) = plotHeight(i);
    yVals(2,i) = plotHeight(i);
    yVals(3,i) = plotHeight(i+1);
    yVals(4,i) = plotHeight(i+1);
    
end % END FOR

patchCoord(:,:,1) = xVals;
patchCoord(:,:,2) = yVals;

end % END FUNCTION

function [ patchCoord ] = GenPatchCoord2( numChans )
% This function will assume 16 channels
if numChans == 17
    numSquare = 16 * 16; % Each electrode will have 16 plots associated with it.
else
    numSquare = numChans * numChans;
end % END IF

patchCoord = zeros(4,4,2);
for i = 1:4
    patchCoord(:,i,1) = (i-1) * 4; % x values
    patchCoord(i,:,2) = 16 - (i-1)*4; % y values
end % END IF

patchCoord = reshape( permute(patchCoord,[2,1,3]), 16, 1, 2);

end % END FUNCTION

% EOF