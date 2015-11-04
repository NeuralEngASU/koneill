function ParametersStruct = parseNEVMarkers(NEV)

% Parses digital data from an unparsed NEV file. The resulting structure
% array contains the timestamps (30 kS/s) and parameters for each trial.
% The last structure of the array contains the database task parameters.
% The length of the array is the number of trials + 1.
%
% Version date: 20121005
% Version changes: Updated to work with video markers
% Author: Tyler Davis

% Parsing NEV digital data
DigitalData = char(NEV.Data.SerialDigitalIO.UnparsedData);
DigitalDataTS = NEV.Data.SerialDigitalIO.TimeStamp;

ParsedData = regexp(DigitalData,'*','split')';
ParsedData = ParsedData(2:end);
ParsedDataTS = double(DigitalDataTS(regexp(DigitalData,'*'))');

DataTypes = regexp(ParsedData,'^(\w+):','tokens');
DataTypes = [DataTypes{:}];
DataTypes = unique([DataTypes{:}]);

MarkersIdxs = false(size(ParsedData,1),1);
for k=1:length(DataTypes)    
    DataTypeIdxs = ~cellfun(@isempty,regexp(ParsedData,['^',DataTypes{k},':']));
    MarkersIdxs = MarkersIdxs | DataTypeIdxs;    
    eval([DataTypes{k},' = ParsedData(DataTypeIdxs);'])
    eval([DataTypes{k},'TS = ParsedDataTS(DataTypeIdxs);'])
end
Markers = ParsedData(~MarkersIdxs);
MarkersTS = ParsedDataTS(~MarkersIdxs);

% Creating a structure containing the data in trial format
TrialParameters = regexprep(TrialParameters,'^TrialParameters:|;#$','');
TrialParameters = regexp(TrialParameters,';','split');
TrialParametersFields = regexp([TrialParameters{:}],'=','split');
TrialParametersFields = reshape([TrialParametersFields{:}],2,[])';
TrialParametersFields = unique(TrialParametersFields(:,1));

TaskParameters = regexprep(cell2mat(TaskParameters),'^TaskParameters:|;#.*$','');
TaskParameters = regexp(TaskParameters,';','split');
TaskParameters = regexp(TaskParameters,'=','split');    
TaskParameters = reshape([TaskParameters{:}],2,[])';
TaskParametersFields = TaskParameters(:,1);

MarkersFields = cellfun(@(x) [x,'TS'],unique(regexprep(Markers,'=|#|[\d\(\)\.]+','')),'uniformoutput',false);
MarkersFields(~cellfun(@isempty,regexp(MarkersFields,'Comments'))) = {'CommentsTS'}; %adding a comments field if it exists
if any(strcmp(MarkersFields,'FramesTS'))
    MarkersFields = ['Frames';MarkersFields];
end
if any(strcmp(MarkersFields,'CommentsTS'))
    MarkersFields = ['Comments';MarkersFields];
end

TrialParametersFields = [MarkersFields;'TrialParametersTS';TrialParametersFields;TaskParametersFields];
TrialParametersFields = regexprep(TrialParametersFields,'/','_');
ParametersStruct = repmat(cell2struct(cell(size(TrialParametersFields,1),1),TrialParametersFields,1),size(TrialParameters,1)+2,1);

% Populating trial parameter fields
for k=1:size(TrialParameters,1)
    TMPParameters = regexp(TrialParameters{k,1},'=','split');    
    TMPParameters = reshape([TMPParameters{:}],2,[])';
    ParametersStruct(k).TrialParametersTS = TrialParametersTS(k);
    for m=1:size(TMPParameters,1)
        ParametersStruct(k).(TMPParameters{m,1}) = TMPParameters{m,2};        
    end    
end

% Populating marker fields
TrialTimes = [ParametersStruct.TrialParametersTS]';
TrialTimes = [[0;TrialTimes(1:end-1)+1],TrialTimes];
TrialTimes = [TrialTimes;[TrialTimes(end,2)+1,TaskParametersTS-1]];
for k=1:size(TrialTimes,1)
    TMPTrialMarkers = regexprep(Markers(MarkersTS>=TrialTimes(k,1) & MarkersTS<=TrialTimes(k,2)),'=|#|[\d\(\)\.]+','');    
    TMPTrialMarkers = cellfun(@(x) [x,'TS'],TMPTrialMarkers,'uniformoutput',false);
    TMPTrialMarkers(~cellfun(@isempty,regexp(TMPTrialMarkers,'Comments'))) = {'CommentsTS'};
    TMPTrialMarkersTS = MarkersTS(MarkersTS>=TrialTimes(k,1) & MarkersTS<=TrialTimes(k,2)); 
    [UniqueTMPTrialMarkers,~,UniqueTMPTrialMarkersIdxs] = unique(TMPTrialMarkers);
    TMPMarkerData = regexp(Markers(MarkersTS>=TrialTimes(k,1) & MarkersTS<=TrialTimes(k,2)),'=','split');
    TMPMarkerData = reshape([TMPMarkerData{:}],2,[])';
    TMPFrames = TMPMarkerData; TMPFrames(cellfun(@isempty,regexp(TMPFrames(:,1),'Frame')),2)={''};
    TMPFrames = regexp(TMPFrames(:,2),'^\d+','match');
    TMPComments = TMPMarkerData; TMPComments(cellfun(@isempty,regexp(TMPComments(:,1),'Comment')),2)={''};
    TMPComments = TMPComments(:,2);
    for m=1:length(UniqueTMPTrialMarkers)
        ParametersStruct(k,1).(regexprep(UniqueTMPTrialMarkers{m},'/','_')) = TMPTrialMarkersTS(UniqueTMPTrialMarkersIdxs==m);
        if strcmp(UniqueTMPTrialMarkers(m),'FramesTS')
            ParametersStruct(k,1).('Frames') = cellfun(@str2num,[TMPFrames{UniqueTMPTrialMarkersIdxs==m}]');
        end
        if strcmp(UniqueTMPTrialMarkers(m),'CommentsTS')
            ParametersStruct(k,1).('Comments') = TMPComments(UniqueTMPTrialMarkersIdxs==m);
        end
    end    
end

% Populating task parameter fields
for m=1:size(TaskParameters,1)
    ParametersStruct(end).(regexprep(TaskParameters{m,1},'/','_')) = TaskParameters{m,2};
end

