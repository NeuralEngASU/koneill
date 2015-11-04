function ParametersStruct = parseNEV_FingerPress_Kevin(NEV)

% Parses digital data from an unparsed NEV file recorded using the
% ClinicalCart ArmMovement.vi task. The resulting structure array contains
% the timestamps (30 kS/s) and parameters for each trial. The last
% structure of the array contains the database task parameters. The length
% of the array is the number of trials + 1.
%
% Version date: 20120514
% Author: Kevin O'Neill
% Original Author: Tyler Davis


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
TrialParametersFields = cellfun(@(x) [x],unique(regexprep(TrialParametersFields,'\.','_')),'uniformoutput',false);

TaskParameters = regexprep(TaskParameters{1},'^TaskParameters:|;#.*$','');
TaskParameters = regexp(TaskParameters,';','split');
TaskParameters = regexp(TaskParameters,'=','split');    
TaskParameters = reshape([TaskParameters{:}],2,[])';
TaskParametersFields = TaskParameters(:,1);
TaskParametersFields = cellfun(@(x) [x],unique(regexprep(TaskParametersFields,'\.','_')),'uniformoutput',false);


MarkersFields = cellfun(@(x) [x,'TS'],unique(regexprep(Markers,'=|#','')),'uniformoutput',false);

markersTest = regexp([MarkersFields{49:end}], '(Frame[0-9]+)', 'Tokens')';
for loop=1:length(markersTest)
    markersTest{loop} = strcat(markersTest{loop}{1}, 'TS');
end % END FOR

markersTest2 = regexp([MarkersFields{1:48}], '([0-9]+)([A-Za-z]+)(TS)', 'Tokens')';
for loop = 1:length(markersTest2)
    markersTest2{loop} = strcat(markersTest2{loop}{2}, markersTest2{loop}{1}, markersTest2{loop}{3});
end % END FOR

MarkersFields = [markersTest2; markersTest];

% MarkersFields = MarkersFields([2,1,4,3]);
% %%%%%%%%%%%%%%%%%
% MarkersFields = cellfun(@(x) [x],unique(regexprep(MarkersFields,'\.','_')),'uniformoutput',false);
% MarkersFields = cellfun(@(x) [x],unique(regexprep(MarkersFields,'\(','_')),'uniformoutput',false);
% MarkersFields = cellfun(@(x) [x],unique(regexprep(MarkersFields,'\)','_')),'uniformoutput',false);
% MarkersFields = cellfun(@(x) ['A', x],MarkersFields, 'uniformoutput',false);
% 
% %%%%%%%%%%%%%%%%%
% MarkersFields = MarkersFields([2,1,4,3]);

TrialParametersFields = [MarkersFields;'TrialParametersTS';TrialParametersFields;TaskParametersFields];
ParametersStruct = repmat(cell2struct(cell(size(TrialParametersFields,1),1),TrialParametersFields,1),size(TrialParameters,1)+1,1);

% for i = 1:length(TrialParameters)
%     TrialParameters{i}{1} = regexprep(TMPParameters{i}{1},'\.','_');
% end % END FOR

% Populating trial parameter fields
for k=1:size(TrialParameters,1)
    TMPParameters = regexp(TrialParameters{k,1},'=','split');   
    for i = 1:length(TMPParameters)
        TMPParameters{i}{1} = regexprep(TMPParameters{i}{1},'\.','_');
    end % END FOR
    TMPParameters = reshape([TMPParameters{:}],2,[])';
    ParametersStruct(k).TrialParametersTS = TrialParametersTS(k);
    for m=1:size(TMPParameters,1)
        ParametersStruct(k).(TMPParameters{m,1}) = TMPParameters{m,2};        
    end    
end

% Populating marker fields
TrialTimes = [ParametersStruct.TrialParametersTS]';
TrialTimes = [[0;TrialTimes(1:end-1)],TrialTimes];
for k=1:size(TrialTimes,1)
    TMPTrialMarkers = regexprep(Markers(MarkersTS>=TrialTimes(k,1) & MarkersTS<TrialTimes(k,2)),'=|#','');
    %%%%%
    TMPTrialMarkers2 = cellfun(@(x) [x],regexp(TMPTrialMarkers,'(Frame[0-9]+)','Tokens'), 'uniformoutput', false);
    TMPTrialMarkers3 = cellfun(@(x) [x],regexp(TMPTrialMarkers,'(^[0-9]+)([A-Za-z]+)','Tokens'), 'uniformoutput', false);
    
    for loop = 1:length(TMPTrialMarkers2)
        
        if isempty(TMPTrialMarkers3{loop})
            TMPTrialMarkers(loop) = TMPTrialMarkers2{loop}{1};
        else
            TMPTrialMarkers{loop} = strcat(TMPTrialMarkers3{loop}{1}{2}, TMPTrialMarkers3{loop}{1}{1});
        end % END IF
        
    end % END FOR
    
    %%%%%
    TMPTrialMarkers = cellfun(@(x) [x,'TS'],TMPTrialMarkers,'uniformoutput',false);
    TMPTrialMarkersTS = MarkersTS(MarkersTS>=TrialTimes(k,1) & MarkersTS<TrialTimes(k,2));            
    for m=1:size(TMPTrialMarkers,1)
%         disp([num2str(m), '/', num2str(size(TMPTrialMarkers,1))])
        ParametersStruct(k,1).(TMPTrialMarkers{m,1}) = TMPTrialMarkersTS(m);
    end
end

TaskParameters(:,1) = cellfun(@(x) [x],regexprep(TaskParameters(:,1)','\.','_'), 'uniformoutput',false);

% Populating task parameter fields
for m=1:size(TaskParameters,1)
    ParametersStruct(end).(TaskParameters{m,1}) = TaskParameters{m,2};
end

