function Header = tank2mat( tank, block, varargin )

Data = struct('Header', [], 'Spikes', [], 'Channels', []);

% defaults
T1       = 0;
T2       = 0;
RANGES   = [];
VERBOSE  = 1;
TYPE     = 1;
SORTNAME = 'TankSort';
SERVER   = 'Local';
NODATA   = false;
CHANNEL  = 0;
STORE    = '';
TTX      = [];

MAXEVENTS = 1e6;
MAXCHANNELS = 1024;

% parse varargin
for i = 1:2:length(varargin)
    eval([upper(varargin{i}) '=varargin{i+1};']);
end

if TYPE == 1, TYPE = 1:5; end
ReadEventsOptions = 'ALL';
if NODATA, ReadEventsOptions = 'NODATA'; end
if ~isscalar(CHANNEL), error('CHANNEL must be a scalar'), end
if CHANNEL < 0, error('CHANNEL must be non-negative'), end
CHANNEL = int32(CHANNEL);

bUseOutsideTTX = ~isempty(TTX);

if ~bUseOutsideTTX
    % create TTankX object
    h = figure('Visible', 'off', 'HandleVisibility', 'off');
    TTX = actxcontrol('TTank.X', 'Parent', h);

    % connect to server
    if TTX.ConnectServer(SERVER, 'TDT2mat') ~= 1
        close(h)
        error(['Problem connecting to server: ' SERVER])
    end

    % open tank
    if TTX.OpenTank(tank, 'R') ~= 1
        TTX.ReleaseServer;
        close(h);
        error(['Problem opening tank: ' tank]);
    end

    % select block
    if TTX.SelectBlock(['~' block]) ~= 1
        block_name = TTX.QueryBlockName(0);
        block_ind = 1;
        while strcmp(block_name, '') == 0
            block_ind = block_ind+1;
            block_name = TTX.QueryBlockName(block_ind);
            if strcmp(block_name, block)
                error(['Block found, but problem selecting it: ' block]);
            end
        end
        error(['Block not found: ' block]);
    end
end

% set info fields
start = TTX.CurBlockStartTime;
stop = TTX.CurBlockStopTime;
total = stop-start;

Data.Header.tankpath = TTX.GetTankItem(tank, 'PT');
Data.Header.blockname = block;
Data.Header.date = TTX.FancyTime(start, 'Y-O-D');
Data.Header.starttime = TTX.FancyTime(start, 'H:M:S');
Data.Header.stoptime = TTX.FancyTime(stop, 'H:M:S');
if stop > 0
    Data.Header.duration = TTX.FancyTime(total, 'H:M:S');
end
Data.Header.streamchannel = CHANNEL;

%data.info.notes = {};

%ind = 1;
%note = TTX.GetNote(ind);
%while ~strcmp(note, '')
%    data.info.notes{ind} = note;
%    ind = ind + 1;
%    note = TTX.GetNote(ind);
%end

if VERBOSE
    fprintf('\nTank Name:\t%s\n', tank);
    fprintf('Tank Path:\t%s\n', Data.Header.tankpath);
    fprintf('Block Name:\t%s\n', Data.Header.blockname);
    fprintf('Start Date:\t%s\n', Data.Header.date);
    fprintf('Start Time:\t%s\n', Data.Header.starttime);
    if stop > 0
        fprintf('Stop Time:\t%s\n', Data.Header.stoptime);
        fprintf('Total Time:\t%s\n', Data.Header.duration);
    else
        fprintf('==Block currently recording==\n');
    end
end

% set global tank server defaults
TTX.SetGlobalV('WavesMemLimit',1e9);
TTX.SetGlobalV('MaxReturn',MAXEVENTS);
TTX.SetGlobalV('T1', T1);
TTX.SetGlobalV('T2', T2);
TTX.SetGlobalV('Channel', CHANNEL);

start = TTX.CurBlockStartTime;
stop = TTX.CurBlockStopTime;
total = stop - start;

fprintf('Start Time: %s\n', TTX.FancyTime(start, 'D/O/Y H:M:U'));
fprintf('Stop Time: %s\n', TTX.FancyTime(stop, 'D/O/Y H:M:U'));
fprintf('Total Time: %s\n\n', TTX.FancyTime(total, 'H:M:U'));

lStores = TTX.GetEventCodes(0);
for i = 1:length(lStores)
    name = TTX.CodeToString(lStores(i));
    fprintf('Store Name:\t%s\n', name)
    
    TTX.GetCodeSpecs(lStores(i));
    type = TTX.EvTypeToString(TTX.EvType);
    fprintf('EvType:\t\t%s\n', type)


    if strcmp(type,'Stream') || bitand(TTX.EvType, 33025) == 33025 % catch RS4 header (33073)
        fprintf('Samp Rate:\t%f\n',TTX.EvSampFreq)
        N = TTX.ReadEventsV(10000, name, 0, 0, 0, 0, 'ALL');
        num_channels = max(TTX.ParseEvInfoV(0, N, 4));
        fprintf('Channels:\t%d\n', num_channels)
        
        for j = 1:num_channels
            TTX.SetGlobalV('Channel', j);
            
            dataTmp = TTX.ReadWavesV(' name ');
            
            h5create('TestDump.h5', ['/C', num2str(j)], [1, length(dataTmp)]);
            h5write ('TestDump.h5', ['/C', num2str(j)], dataTmp);
                    
        end % END FOR
        
        clear('dataTmp');
        
        Header.fs = TTX.EvSampFreq;

    end
    disp(' ')

end % END FUNCTION

% EOF