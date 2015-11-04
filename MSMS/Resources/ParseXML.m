%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:    ParseXML
% Desc:     Parses the "simulationSetup.xml" for MSMS
% Author:   Kevin O'Neill (Greger Lab)
% Date:     July 1, 2012
% 
% Use:      myStruct = ParseXML();   -- Will prompt user with a uigetfile()
%           mystruct = ParseXML('myFileName.xml')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% E:\MSMS\MSMS 2.1\Tutorials\FingerTaskDemo\Model\Prostheses
% E:\MSMS\MSMS 2.1\Tutorials\FingerTaskDemo\Setup

function SimInfo = ParseXML( folderName )
% ParseXML: Parses simulationSetup.xml
%   Parses the XML file into a structure with type-correct values

if nargin == 1
    try
        setupXML  = fileread([folderName, '\Setup\simulationSetup.xml']);
        offsetXML = fileread([folderName, '\Model\Prostheses\prosthesis.xml']);

        % Opens the JoinFile.ini which contains more natural names for joints
        jointFile = fileread('JointFile.ini');
    catch
        error('Failed to read XML file: %s.\nFile may not exist or is not in MATLAB''s path.',fileName);
    end % END TRY
else
    % Prompts the user to select an 'scenario' for MSMS.
    folderName = uigetdir();
    % Opens the files and splits it into a cell based on the newline character
    try
%         setupXML  = fileread([folderName, '\Setup\simulationSetup.xml']);
%         offsetXML = fileread([folderName, '\Model\Prostheses\prostheses.xml']);
        setupXML  = fileread(['simulationSetup.xml']);
        offsetXML = fileread(['prosthesis.xml']);
        % Opens the JoinFile.ini which contains more natural names for joints
        jointFile = fileread('JointFile.ini'); 
    catch
        error('Failed to read files: %s.\nFiles may not exist, files may not be in MATLAB''s path, or wrong folder was selected.',folderName);
    end % END TRY
end

setupXMLSplit = regexp(setupXML, '\n', 'split')';

% Matches the XML type and encoding options
matchXMLInfo = regexp(setupXMLSplit{1},'"(.*?)"','tokens');

% Stores the XML version, encoding, and standalone into SimInfo
SimInfo.XML.version    = str2double(matchXMLInfo{1});
SimInfo.XML.encoding   = char(matchXMLInfo{2});
SimInfo.XML.standalone = char(matchXMLInfo{3});

% Matches the Setupinfo type and encoding options
matchSetupInfo = regexp(setupXMLSplit{2},'"(.*?)"','tokens');

% Stores the XML version, encoding, and standalone into SimInfo
SimInfo.SetupInfo.creationDate = char(matchSetupInfo{1});
SimInfo.SetupInfo.setupName    = char(matchSetupInfo{2});
SimInfo.SetupInfo.version      = str2double(matchSetupInfo{3});
SimInfo.SetupInfo.xmlns        = char(matchSetupInfo{4});

% Expresions for tabs and components
tabsExpr{1} = '<general>.*?</general>';
tabsExpr{2} = '<setup>.*?</setup>';
tabsExpr{3} = '<solver>.*?</solver>';
tabsExpr{4} = '<dynamicEngine>.*?</dynamicEngine>'; 
tabsExpr{5} = '<outputData>.*?</outputData>';

componentExpr = '<component>.*?</component>';

% Cycle through the text file looking for all of the tabs and place the
% tabs into a cell
tabsCell = cell(length(tabsExpr));
for k = 1 : length(tabsExpr)
    tabsCell{k} = regexp(setupXML,tabsExpr{k},'match')';
end % END FOR

% Look through the text file to find all of the components
componentCell = regexp(setupXML,componentExpr,'match')';

% Parse all of the tabs and components and assign to SimulationStruct
GeneralStruct       = ParseGeneral(tabsCell{1});
ComponentsStruct    = ParseComponent(componentCell, jointFile, offsetXML);
SolverStruct        = ParseSolver(tabsCell{3});
DynamicEngineStruct = ParseDynamicEngine(tabsCell{4});
OutoutDataStruct    = ParseOutputData(tabsCell{5});

% Assigns the parsed structs to the returned struct (This workflow was
% chosen for debugging purposes)
SimInfo.General       = GeneralStruct;
SimInfo.Components    = ComponentsStruct;
SimInfo.Solver        = SolverStruct;
SimInfo.DynamicEngine = DynamicEngineStruct;
SimInfo.OutputData    = OutoutDataStruct;

SimInfo.Components = orderfields(SimInfo.Components);

end % END FUNCTION


%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function GeneralStruct = ParseGeneral(generalTab)
% Parses the general tab of the XML file

% Converts the cell into a string
generalTab = char(generalTab);

% Expressions to take out exportFilename and gravity information
generalExpr{1} = '<exportFilename>(.+?)</exportFilename>';
generalExpr{2} = '<(ns[0-9]+):(.).+?"(.+?)">(.+?)<';

% Extract the information from string
exportFilename = regexp(generalTab,generalExpr{1},'tokens');
gravity        = regexp(generalTab,generalExpr{2},'tokens');

% Assign exportFilename tab to structure (assumed to be logical)
GeneralStruct.ExportFilename = logical(str2double(exportFilename{1}));

% Loops through all given directions and assign gravity information
for k = 1:length(gravity)
    eval(['GeneralStruct.Gravity.', gravity{k}{2}, '.ns = ''', gravity{k}{1}, ''';']);
    eval(['GeneralStruct.Gravity.', gravity{k}{2}, '.Value = str2double(''', gravity{k}{4}, ''');']);
    eval(['GeneralStruct.Gravity.', gravity{k}{2}, '.xmlns = ''', gravity{k}{3}, ''';']);
end % END FOR
end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function ComponentStruct = ParseComponent(componentCell, jointFile, offsetXML)
% Parses the component(s) tab of the XML file

% Expressions to take out subcomponent information
componentExpr{1} = '<jointName>(.+?)</jointName>';
componentExpr{2} = '<seqNum>(.+?)</seqNum>';
componentExpr{3} = '<simulationType>(.+?)</simulationType>';
componentExpr{4} = '<name>(.+?)</name>';
componentExpr{5} = '<initialPositionMOrRad>(.+?)</initialPositionMOrRad>';
componentExpr{6} = '<velocityMeterPerSecOrRadPerSec>(.+?)</velocityMeterPerSecOrRadPerSec>';
componentExpr{7} = '<(ns[0-9]+):(.+?)\s.+?"(.+?)">(.+?)<';

% Loops over each component
for k = 1:length(componentCell)
    % Extract the information from string
    jointName = regexp(componentCell{k},componentExpr{1},'tokens');
    seqNum    = regexp(componentCell{k},componentExpr{2},'tokens');
    simType   = regexp(componentCell{k},componentExpr{3},'tokens');
    nameDOF   = regexp(componentCell{k},componentExpr{4},'tokens');
    initPos   = regexp(componentCell{k},componentExpr{5},'tokens');
    velocity  = regexp(componentCell{k},componentExpr{6},'tokens');
    jointID   = regexp(componentCell{k},componentExpr{7},'tokens');
    
    jointStructName = ParseJointName(jointName, jointFile);
    
    % Loops through all DOFs that a component has and assigns them to a struct
    for j = 1:length(nameDOF)
        eval(['ComponentStruct.', jointStructName, '.DOF', num2str(j),'.NameDOF = ''', char(nameDOF{j}), ''';']);
        eval(['ComponentStruct.', jointStructName, '.DOF', num2str(j),'.InitPos = ', char(initPos{j}), ';']);
        eval(['ComponentStruct.', jointStructName, '.DOF', num2str(j),'.Velocity = ', char(velocity{j}), ';']);
    end % END FOR
    
    % Assigns the number of DOFs to nDOF
    if ~isempty(j)
        eval(['ComponentStruct.', jointStructName, '.nDOF = ', num2str(j),';']);
    end
    
    % If the characterID is one character long, add a white space as the
    % preceeding character
    if length(jointID{1}{4}) == 1
        jointID{1}{4} = [' ', jointID{1}{4}];
    end
    
    % Assigns all componentType information for jointId
    eval(['ComponentStruct.', jointStructName, '.', jointID{1}{2}, '.ns= ''',  jointID{1}{1}, ''';']);
    eval(['ComponentStruct.', jointStructName, '.', jointID{1}{2}, '.charID= ''',  jointID{1}{4}, ''';']);
    eval(['ComponentStruct.', jointStructName, '.', jointID{1}{2}, '.xmlns= ''',  jointID{1}{3}, ''';']);
    
    % Assigns all componentNumber information for jointId
    eval(['ComponentStruct.', jointStructName, '.', jointID{2}{2}, '.ns= ''',  jointID{2}{1}, ''';']);
    eval(['ComponentStruct.', jointStructName, '.', jointID{2}{2}, '.xmlns= ''',  jointID{2}{3}, ''';']);
    eval(['ComponentStruct.', jointStructName, '.', char(jointID{2}{2}), '.numID= int16(',  char(jointID{2}{4}), ');']);

    % Assigns the jointName to the struct
    eval(['ComponentStruct.', jointStructName, '.Name = ''',  char(jointName{1}), ''';']);

    % If the sequence number exists, assigns it to the struct
    if ~isempty(seqNum)
        eval(['ComponentStruct.', jointStructName, '.SeqNum = ',  char(seqNum{1}), ';']);
    end
    
    % Assigns the simulation type to the struct
    eval(['ComponentStruct.', jointStructName, '.SimType = ''',  char(simType{1}), ''';']);       
end % END FOR

offsetExpr{1} = '<component xsi:type="[A-Za-z]+Joint" name="(.+?)" xmlns:xsi=".+?">(.+?)</component>';
offsetExpr{2} = '<inboardJointCenterM>\s+<x>(.+?)</x>\s+<y>(.+?)</y>\s+<z>(.+?)</z>\s+</inboardJointCenterM>';
offsetExpr{3} = '<axisM>\s+<x>(.+?)</x>\s+<y>(.+?)</y>\s+<z>(.+?)</z>\s+</axisM>';

offsetComponent = regexp(offsetXML,offsetExpr{1},'Tokens');

for k = 1:length(offsetComponent)
    
    offsetABC = regexp(offsetComponent{k}{2}, offsetExpr{2}, 'Tokens');
    offsetUVW = regexp(offsetComponent{k}{2}, offsetExpr{3}, 'Tokens');
    
    jointStructName = ParseJointName(offsetComponent{k}{1}, jointFile);
    
    offsetABCStr = ['[', offsetABC{1}{1}, ',', offsetABC{1}{2}, ',', offsetABC{1}{3}, ']'];
    
    if ~isempty(offsetUVW)
        offsetUVWStr = ['[', offsetUVW{1}{1}, ',', offsetUVW{1}{2}, ',', offsetUVW{1}{3}, ']'];
    else
        offsetUVWStr = 'NaN';
    end % END IF
    
    eval(['ComponentStruct.', jointStructName, '.abc = ',  offsetABCStr,';']);
    eval(['ComponentStruct.', jointStructName, '.uvw = ',  offsetUVWStr,';']);

end % END FOR
end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function jointStructName = ParseJointName(jointName, jointFile)

% Converts the cell into a string
if iscell(jointName)
    jointName = char(jointName{1});
end

jointNameExpr = [jointName, ':([A-Za-z0-9]+)'];

jointStructName = regexp(jointFile, jointNameExpr,'tokens');
jointStructName = char(jointStructName{1});

end

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function SolverStruct = ParseSolver(solverTab)
% Parses the solver tab of the XML file

% Converts the cell into a string
solverTab = char(solverTab);

% Expressions to take out solver information
solverExpr{1}  = '<startTimeSec>(.+?)</startTimeSec>';
solverExpr{2}  = '<stopTimeSec>(.+?)</stopTimeSec>';
solverExpr{3}  = '<solverName>(.+?)</solverName>';
solverExpr{4}  = '<integrationType>(.+?)</integrationType>';
solverExpr{5}  = '<fixed>(.+?)</fixed>';
solverExpr{6}  = '<min>(.+?)</min>';
solverExpr{7}  = '<max>(.+?)</max>';
solverExpr{8}  = '<initial>(.+?)</initial>';
solverExpr{9}  = '<absolute>(.+?)</absolute>';
solverExpr{10} = '<relative>(.+?)</relative>';

% Extracts solver information
startTime  = regexp(solverTab,solverExpr{1}, 'tokens');
stopTime   = regexp(solverTab,solverExpr{2}, 'tokens');
solverName = regexp(solverTab,solverExpr{3}, 'tokens');
intType    = regexp(solverTab,solverExpr{4}, 'tokens');
fixedSec   = regexp(solverTab,solverExpr{5}, 'tokens');
minSec     = regexp(solverTab,solverExpr{6}, 'tokens');
maxSec     = regexp(solverTab,solverExpr{7}, 'tokens');
initialSec = regexp(solverTab,solverExpr{8}, 'tokens');
absTol     = regexp(solverTab,solverExpr{9}, 'tokens');
relTol     = regexp(solverTab,solverExpr{10},'tokens');

% Assigns the solver information to the struct
SolverStruct.SimulationTime.StartTime          = str2double(char(startTime{1}));
SolverStruct.SimulationTime.StopTime           = str2double(char(stopTime{1}));
SolverStruct.SolverOptions.SolverName          = char(solverName{1});
SolverStruct.SolverOptions.IntegrationType     = char(intType{1});
SolverStruct.SolverOptions.StepSizeSec.Fixed   = str2double(char(fixedSec{1}));
SolverStruct.SolverOptions.StepSizeSec.Min     = str2double(char(minSec{1}));
SolverStruct.SolverOptions.StepSizeSec.Max     = str2double(char(maxSec{1}));
SolverStruct.SolverOptions.StepSizeSec.Initial = str2double(char(initialSec{1}));
SolverStruct.SolverOptions.Tolerance.Absolute  = str2double(char(absTol{1}));
SolverStruct.SolverOptions.Tolerance.Relative  = str2double(char(relTol{1}));

end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function DynamicEngineStruct = ParseDynamicEngine(dynamicEngineTab)
% Parses the dynamicEngine tab of the XML file

% Converts the cell into a string
dynamicEngineTab = char(dynamicEngineTab);

% Extracts the name of the DynamicEngine
dynamicEngineExpr = '<name>(.+?)</name>';
engineName = regexp(dynamicEngineTab,dynamicEngineExpr,'tokens');

% Assigns the name of the engine to the struct
DynamicEngineStruct.EngineName = char(engineName{1});

end % END FUNCTION

%%%%%%%%%%%%%%%% - SUBFUNCTION - %%%%%%%%%%%%%%%%
function OutputDataStruct = ParseOutputData(outputDataTab)
% Parses the outputData tab of the XML file

% Converts the cell into a string
outputDataTab = char(outputDataTab);

% Expressions to take out outputData information
outputDataExpr{1} = '<filename>(.+?)</filename>';
outputDataExpr{2} = '<samplingTimeSec>(.+?)</samplingTimeSec>';
outputDataExpr{3} = '<dataStreamingMethod>(.+?)</dataStreamingMethod>';
outputDataExpr{4} = '<animationDataMethod>(.+?)</animationDataMethod>';

% Extracts outputData information
fileName  = regexp(outputDataTab,outputDataExpr{1},'tokens');
samplingTime   = regexp(outputDataTab,outputDataExpr{2},'tokens');
dataStream = regexp(outputDataTab,outputDataExpr{3},'tokens');
aniData    = regexp(outputDataTab,outputDataExpr{4},'tokens');

% Assigns the outputData information to the struct
OutputDataStruct.FileStorage.FileName              = char(fileName{1});
OutputDataStruct.FileStorage.SamplingTimeSec       = str2double(char(samplingTime{1}));
OutputDataStruct.LiveAnimation.DataStreamingMethod = char(dataStream{1});
OutputDataStruct.LiveAnimation.AnimationDataMethod = char(aniData{1});
OutputDataStruct.LiveAnimation.SamplingTimeSec     = str2double(char(samplingTime{2}));

end % END FUNCTION
% EOF