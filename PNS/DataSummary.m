function varargout = DataSummary(varargin)
% DATASUMMARY MATLAB code for DataSummary.fig
%      DATASUMMARY, by itself, creates a new DATASUMMARY or raises the existing
%      singleton*.
%
%      H = DATASUMMARY returns the handle to a new DATASUMMARY or the handle to
%      the existing singleton*.
%
%      DATASUMMARY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATASUMMARY.M with the given input arguments.
%
%      DATASUMMARY('Property','Value',...) creates a new DATASUMMARY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DataSummary_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DataSummary_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DataSummary

% Last Modified by GUIDE v2.5 06-Jul-2012 10:22:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DataSummary_OpeningFcn, ...
                   'gui_OutputFcn',  @DataSummary_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before DataSummary is made visible.
function DataSummary_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DataSummary (see VARARGIN)

% Choose default command line output for DataSummary
handles.output = hObject;

%%% Load Panel %%%
% Initial Values
handles.Flag.neuralDataCheck  = false;
handles.Flag.voltageDataCheck = false;
handles.Flag.markerDataCheck  = false;
handles.Flag.markerOnSpect    = false;
handles.Flag.markerOnVoltage  = false;

handles.Inputs.neuralDataChannelRange = '1:64';
handles.Inputs.voltageDataChannelRange = '142,143,144';

handles.Data.neuralData  = [];
handles.Data.voltageData = [];
handles.Data.markerData  = [];

% GUI setup
set(handles.NeuralDataCheck,  'Value', handles.Flag.neuralDataCheck);
set(handles.VoltageDataCheck, 'Value', handles.Flag.voltageDataCheck);
set(handles.MarkerDataCheck,  'Value', handles.Flag.markerDataCheck);

set(handles.NeuralDataChannelRange, 'String', handles.Inputs.neuralDataChannelRange);
set(handles.NeuralDataChannelRange, 'Enable', 'Off');

set(handles.VoltageDataChannelRange, 'String', handles.Inputs.voltageDataChannelRange);
set(handles.VoltageDataChannelRange, 'Enable', 'Off');

set(handles.NeuralDataLoad,  'Enable', 'Off')
set(handles.VoltageDataLoad, 'Enable', 'Off')
set(handles.MarkerDataLoad,  'Enable', 'Off')
set(handles.MarkerOnSpect,   'Enable', 'Off')
set(handles.MarkerOnVoltage, 'Enable', 'Off')

%%% Options Panel %%%
handles.Flag.removeArtifacts    = false;
handles.Flag.interquartileRange = false;
handles.Flag.averageVoltage     = false;

handles.Inputs.avgVoltSecRangeMin = 0.0;
handles.Inputs.avgVoltSecRangeMax = 1.0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DataSummary wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DataSummary_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - Load Panel - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in NeuralDataCheck.
function NeuralDataCheck_Callback(hObject, eventdata, handles)

% Gets the current value from the checkbox and assigns it to neuralDataCheck
handles.Flag.neuralDataCheck = get(hObject, 'Value');

% If the checkbox is checked, enable the pushbutton
if handles.Flag.neuralDataCheck
    set(handles.NeuralDataLoad,         'Enable', 'On')
    set(handles.NeuralDataChannelRange, 'Enable', 'On');
else
    set(handles.NeuralDataLoad,         'Enable', 'Off')
    set(handles.NeuralDataChannelRange, 'Enable', 'Off');
end
guidata(hObject, handles);

% --- Executes on button press in NeuralDataLoad.
function NeuralDataLoad_Callback(hObject, eventdata, handles)

[fileName filePath] = uigetfile('.ns4mat');

handles.Data.rawNeural = fileread([filePath fileName]);
set(handles.NeuralDataStatic, 'String', fileName);

guidata(hObject, handles);

function NeuralDataChannelRange_Callback(hObject, eventdata, handles)

handles.Inputs.neuralDataChannelRange = get(hObject,'String');

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function NeuralDataChannelRange_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in VoltageDataCheck.
function VoltageDataCheck_Callback(hObject, eventdata, handles)

% Gets the current value from the checkbox and assigns it to voltageDataCheck
handles.Flag.voltageDataCheck = get(hObject, 'Value');

% If the checkbox is checked, enable the pushbutton
if handles.Flag.voltageDataCheck
    set(handles.VoltageDataLoad,  'Enable', 'On')
else
    set(handles.VoltageDataLoad,  'Enable', 'Off')
end
guidata(hObject, handles);

% --- Executes on button press in VoltageDataLoad.
function VoltageDataLoad_Callback(hObject, eventdata, handles)

[fileName filePath] = uigetfile('.ns3mat');

handles.Data.rawVoltage = fileread([filePath fileName]);
set(handles.VoltageDataStatic, 'String', fileName);

guidata(hObject, handles);

function VoltageDataChannelRange_Callback(hObject, eventdata, handles)

handles.Inputs.voltageDataChannelRange = get(hObject,'String');

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function VoltageDataChannelRange_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MarkerDataCheck.
function MarkerDataCheck_Callback(hObject, eventdata, handles)

% Gets the current value from the checkbox and assigns it to voltageDataCheck
handles.Flag.markerDataCheck = get(hObject, 'Value');

% If the checkbox is checked, enable the pushbutton
if handles.Flag.markerDataCheck
    set(handles.MarkerDataLoad,  'Enable', 'On')
    set(handles.MarkerOnSpect,   'Enable', 'On')
    set(handles.MarkerOnVoltage, 'Enable', 'On')
else
    set(handles.MarkerDataLoad,  'Enable', 'Off')
    set(handles.MarkerOnSpect,   'Enable', 'Off')
    set(handles.MarkerOnVoltage, 'Enable', 'Off')
end
guidata(hObject, handles);

% --- Executes on button press in MarkerDataLoad.
function MarkerDataLoad_Callback(hObject, eventdata, handles)

[fileName filePath] = uigetfile('.mat');

handles.Data.rawMarker = fileread([filePath fileName]);
set(handles.MarkerDataStatic, 'String', fileName);

guidata(hObject, handles);

% --- Executes on button press in MarkerOnSpect.
function MarkerOnSpect_Callback(hObject, eventdata, handles)

handles.Flag.markerOnSpect = get(hObject, 'Value');

guidata(hObject, handles);

% --- Executes on button press in MarkerOnVoltage.
function MarkerOnVoltage_Callback(hObject, eventdata, handles)

handles.Flag.markerOnVoltage = get(hObject, 'Value');

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% - Options Panel - %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in RemoveArtifacts.
function RemoveArtifacts_Callback(hObject, eventdata, handles)

% --- Executes on button press in InterquartileRange.
function InterquartileRange_Callback(hObject, eventdata, handles)

% --- Executes on button press in AverageVoltage.
function AverageVoltage_Callback(hObject, eventdata, handles)

function AvgVoltRangeMin_Callback(hObject, eventdata, handles)
handles.Inputs.avgVoltSecRangeMin = str2double(get(hObject, 'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function AvgVoltRangeMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function AvgVoltRangeMax_Callback(hObject, eventdata, handles)
handles.Inputs.avgVoltSecRangeMax = str2double(get(hObject, 'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function AvgVoltRangeMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - FUNCTIONS - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% EOF



