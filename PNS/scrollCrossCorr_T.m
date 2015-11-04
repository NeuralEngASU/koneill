function varargout = scrollCrossCorr_T(varargin)
% SCROLLCROSSCORR MATLAB code for scrollCrossCorr.fig
%      SCROLLCROSSCORR, by itself, creates a new SCROLLCROSSCORR or raises the existing
%      singleton*.
%
%      H = SCROLLCROSSCORR returns the handle to a new SCROLLCROSSCORR or the handle to
%      the existing singleton*.
%
%      SCROLLCROSSCORR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCROLLCROSSCORR.M with the given input arguments.
%
%      SCROLLCROSSCORR('Property','Value',...) creates a new SCROLLCROSSCORR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before scrollCrossCorr_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to scrollCrossCorr_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help scrollCrossCorr

% Last Modified by GUIDE v2.5 01-May-2013 05:43:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @scrollCrossCorr_T_OpeningFcn, ...
                   'gui_OutputFcn',  @scrollCrossCorr_T_OutputFcn, ...
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


% --- Executes just before scrollCrossCorr is made visible.
function scrollCrossCorr_T_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to scrollCrossCorr (see VARARGIN)

% filename = 'D:\Tyler\data\PNS\P201301\20130415-120054\20130415-120054-001';
filename = 'D:\Tyler\data\PNS\P201301\20130411-112540\20130411-112540-001';
nsx2mat([filename,'.ns5']);
NEV = openNEV([filename,'.nev']);
%load([filename,'.mat'])
load([filename,'.ns5mat'],'-mat','Header')
DOFIdxs = 137:144;
handles.DOF = {'Thumb','Index','Middle','Ring','Little','IndexIntrinsic','RingIntrinsic','LittleIntrinsic'};
if exist([filename,'_TDataDS.mat'],'file')
    load([filename,'_TDataDS.mat'])
else
    Hd = design(fdesign.lowpass('N,F3dB',4,15,30000),'butter');
    TDataDS = zeros(8,ceil(Header.ChannelLengthSamples/500));
    for k=1:length(DOFIdxs)
        TData = readNSxMat([filename,'.ns5mat'],DOFIdxs(k));
        TData = filter(Hd,TData.data);
        TDataDS(k,:) = TData(1:500:end); %downsample to 60 S/s
    end
    clear('TData')        
    save([filename,'_TDataDS.mat'],'TDataDS')
end
handles.TDataDS = TDataDS/1000; %normalize to 1

if exist([filename,'_RateDS.mat'],'file')
    load([filename,'_RateDS.mat'])
else
    boxwin = 0.3; %boxcar window in sec
    TSCell = cell(96,1);
    WFCell = cell(96,1);
    RateDS = zeros(96,ceil(Header.ChannelLengthSamples/500)); %total length for 60 S/sec
    Hd = design(fdesign.lowpass('N,F3dB',4,15,300),'butter');
    for k=1:96
        clc, disp(k)
        TS = double(NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Electrode==k)); 
        TSCell(k) = {ceil(TS./500)}; %timestamps at 60 S/sec
        WFCell(k) = {double(NEV.Data.Spikes.Waveform(:,NEV.Data.Spikes.Electrode==k))};
        RateB = zeros(1,ceil(Header.ChannelLengthSamples/100)); %total length for 300 S/sec
        if ~isempty(TS)
            RateB(ceil(TS./100)) = 1; %timestamps at 300 S/sec
        end
        %     RateB = conv2(RateB,ones(1,boxwin*60)./boxwin,'same');
        RateB = conv2(RateB,gausswin(boxwin*300)'./(boxwin/2),'same'); %std = N/5 (i.e. std = 90/5 = 18samples = 60ms)
        RateB = filter(Hd,RateB); %filter at 15Hz
        RateDS(k,:) = RateB(1:5:end); %downsample from 300 S/sec to 60 S/sec
    end
    clear('RateB')        
    save([filename,'_RateDS.mat'],'TSCell','WFCell','RateDS')
end
handles.TSCell = TSCell;
handles.WFCell = WFCell; %unit waveforms
handles.RateDS = RateDS/50; %normalize so 50Hz is 1

handles.e = 2; %current electrode
handles.plotwin = 7200; %for 60 S/s
handles.kmax = floor((size(handles.TDataDS,2)-handles.plotwin));
handles.k = floor(get(handles.slider1,'Value')*handles.kmax);
handles.plotidxs = (1:handles.plotwin)+handles.k;

handles.mapidxs = c2e(1:96,'pns');
[~,handles.currdof] = max(abs(mean(handles.TDataDS(:,handles.plotidxs),2)));
handles = calcCorr(handles);
handles.CCH = imagesc(handles.CCMat,'parent',handles.axes1);
handles.CH = imagesc(handles.CMat,'parent',handles.axes2);
set(handles.axes1,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
set(handles.axes2,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
set(get(handles.axes1,'title'),'string','Correlation Coefficient')
set(get(handles.axes2,'title'),'string','Covariance')

handles.pH = plot(handles.axes3,handles.plotidxs,handles.RateDS(e2c(handles.e,'pns'),handles.plotidxs),'k',handles.plotidxs,(handles.TDataDS(:,handles.plotidxs))','linewidth',3);
set(handles.axes3,'xlim',[handles.plotidxs(1),handles.plotidxs(end)],'xtick',handles.plotidxs(rem(handles.plotidxs,3600)==0),'xticklabel',round(handles.plotidxs(rem(handles.plotidxs,3600)==0)/3600),'ylim',[-1,1],'ytick',0:0.5:1,'yticklabel',0:25:50)
set(get(handles.axes3,'title'),'string',[handles.DOF{handles.currdof},', e',num2str(handles.e),', c',num2str(e2c(handles.e,'pns'))])
set(get(handles.axes3,'ylabel'),'string','Rate (Hz)')

wfidxs = handles.TSCell{e2c(handles.e,'pns')}>=handles.plotidxs(1) & handles.TSCell{e2c(handles.e,'pns')}<=handles.plotidxs(end);
wfs = handles.WFCell{e2c(handles.e,'pns')}(:,wfidxs); 
if size(wfs,2)>50
    WFS = wfs(:,fix(linspace(1,size(wfs,2),50)));
else
    WFS = nan(48,50);
    WFS(:,1:size(wfs,2)) = wfs;
end
handles.WFH = plot(handles.axes4,1:48,WFS,'b',1:48,mean(wfs,2),'k');
set(handles.WFH(end),'linewidth',3)
set(handles.axes4,'xlim',[1,48],'ylim',[-400,400],'xtick',[],'xticklabel',[],'ytick',-400:400:400,'yticklabel',-100:100:100)
set(get(handles.axes4,'title'),'string',['e',num2str(handles.e),', c',num2str(e2c(handles.e,'pns'))])

set(handles.slider1,'sliderstep',[0.001,0.01])
handles.slider1_listener = addlistener(handles.slider1,'ContinuousValueChange',@(hObject,eventdata)scrollCrossCorr_T('slider1_Callback',hObject,eventdata,guidata(hObject)));

% Choose default command line output for scrollCrossCorr
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes scrollCrossCorr wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = scrollCrossCorr_T_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.k = floor(get(handles.slider1,'Value')*handles.kmax);
handles.plotidxs = (1:handles.plotwin)+handles.k;
[~,handles.currdof] = max(abs(mean(handles.TDataDS(:,handles.plotidxs),2)));
handles = calcCorr(handles);
plotData(handles);
guidata(hObject, handles);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


function plotData(handles)

set(handles.CCH,'cdata',handles.CCMat)
set(handles.CH,'cdata',handles.CMat)
set(handles.pH(1),'xdata',handles.plotidxs,'ydata',handles.RateDS(e2c(handles.e,'pns'),handles.plotidxs))
for m=1:8
    set(handles.pH(m+1),'xdata',handles.plotidxs,'ydata',handles.TDataDS(m,handles.plotidxs))
end
set(get(handles.axes3,'title'),'string',[handles.DOF{handles.currdof},', e',num2str(handles.e),', c',num2str(e2c(handles.e,'pns'))])
set(handles.axes3,'xlim',[handles.plotidxs(1),handles.plotidxs(end)],'xtick',handles.plotidxs(rem(handles.plotidxs,3600)==0),'xticklabel',round(handles.plotidxs(rem(handles.plotidxs,3600)==0)/3600),'ylim',[-1,1],'ytick',0:0.5:1,'yticklabel',0:25:50)

wfidxs = handles.TSCell{e2c(handles.e,'pns')}>=handles.plotidxs(1) & handles.TSCell{e2c(handles.e,'pns')}<=handles.plotidxs(end);
wfs = handles.WFCell{e2c(handles.e,'pns')}(:,wfidxs); 
if size(wfs,2)>50
    WFS = wfs(:,fix(linspace(1,size(wfs,2),50)));
else
    WFS = nan(48,50);
    WFS(:,1:size(wfs,2)) = wfs;
end
for m=1:50
    set(handles.WFH(m),'ydata',WFS(:,m))
end
set(handles.WFH(51),'ydata',mean(wfs,2))
set(get(handles.axes4,'title'),'string',['e',num2str(handles.e),', c',num2str(e2c(handles.e,'pns'))])


function handles = calcCorr(handles)

plotidxs = handles.plotidxs;
mTDataDS = handles.TDataDS(handles.currdof,plotidxs) - repmat(mean(handles.TDataDS(handles.currdof,plotidxs),2),1,length(plotidxs));
sTDataDS = std(handles.TDataDS(handles.currdof,plotidxs),0,2);
mRateDS = handles.RateDS(:,plotidxs) - repmat(mean(handles.RateDS(:,plotidxs),2),1,length(plotidxs));
sRateDS = std(handles.RateDS(:,plotidxs),0,2);

C = abs((mRateDS*mTDataDS')/(length(plotidxs)-1)); %covariance
CMat = nan(100,1);
CMat(handles.mapidxs) = C(:,1);
handles.CMat = rot90(reshape(CMat,10,10));

CC = C./(sRateDS*sTDataDS'); %pearson's correlation coefficient
CC(isnan(CC))=0; CC(isinf(CC))=0;
CCMat = nan(100,1);
CCMat(handles.mapidxs) = CC(:,1);
handles.CCMat = rot90(reshape(CCMat,10,10));


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% cp1 = fix(get(handles.axes1,'currentpoint')+0.5); x1 = cp1(1,1); y1 = 10-cp1(1,2);
% cp2 = fix(get(handles.axes2,'currentpoint')+0.5); x2 = cp2(1,1); y2 = 10-cp2(1,2);
% if x1>=1 && x1<=10 && y1>=0 && y1<=9
%     disp(1)
% elseif x2>=1 && x2<=10 && y2>=0 && y2<=9
%     disp(2)
% end


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cp1 = fix(get(handles.axes1,'currentpoint')+0.5); x1 = cp1(1,1); y1 = 10-cp1(1,2);
cp2 = fix(get(handles.axes2,'currentpoint')+0.5); x2 = cp2(1,1); y2 = 10-cp2(1,2);
if x1>=1 && x1<=10 && y1>=0 && y1<=9
    e = x1+(y1*10);
elseif x2>=1 && x2<=10 && y2>=0 && y2<=9
    e = x2+(y2*10);
end
if ~any(e==[1,11,81,91])
    handles.e = e;
end
plotData(handles);
guidata(hObject, handles);
