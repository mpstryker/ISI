function varargout = ImanMapAns(varargin)
%IMANMAPANS M-file for ImanMapAns.fig
%      IMANMAPANS, by itself, creates a new IMANMAPANS or raises the existing
%      singleton*.
%
%      H = IMANMAPANS returns the handle to a new IMANMAPANS or the handle to
%      the existing singleton*.
%
%      IMANMAPANS('Property','Value',...) creates a new IMANMAPANS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ImanMapAns_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      IMANMAPANS('CALLBACK') and IMANMAPANS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in IMANMAPANS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ImanMapAns

% Last Modified by GUIDE v2.5 21-Jun-2006 18:05:50

%% On Jun 21, 2006,
%% Now using the Guide-callbacks instead of manually -written ones (IMAcallback.m)
%% All old files saved in the folder of old-Jun2006

%% Designed by Jianhua Cang, last modified on Jun-21-2006

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImanMapAns_OpeningFcn, ...
                   'gui_OutputFcn',  @ImanMapAns_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before ImanMapAns is made visible.
function ImanMapAns_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for ImanMapAns
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ImanMapAns wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%% To initialize, is this necessary? %% JC added on Jun-21-2006
global IMAPhase; 
global mm_per_pixel deg_per_stimcycle filter_size;
global TemplateMethod TemplateValue;
global MapAnsColorMap;

deg_per_stimcycle = str2num(get(handles.spatial_freq_edit, 'string'));
mm_per_pixel = str2num(get(handles.um_per_pixel_edit, 'string'))/1000;
filter_size = str2num(get(handles.filter_size_edit, 'string'));
TemplateMethod = get(handles.TemplateMethod, 'value');
TemplateValue = str2num(get(handles.TemplateValue, 'string'));

% load the colortable
S = load('MapAnsColorWorkspace'); 
MapAnsColorMap = S.MapAnsColor;

    % see the following notes from E:\Retinotopy\MapAnsColor.m
    % %%%%%% using MapAns Colortable
    % %%%% Jianhua Cang, 07-27-05
    % 
    % %% import data from E:\Retinotopy\colortable2.txt 
    % %%%%%%%%%notes from MPS:
    % % Colortable2.txt is the one with 360 lines showing the rgb values for 
    % % each degree of phase at the specified lightness (all 0.5) and the 
    % % specified saturation (all 1.0).
    % % 
    % % So you can get the rgb entries for 256 steps by finding the N*256/360 
    % % entry in eqch table.
    % 
    % MapAnsColorTable = colortable2(:, 4:6)/255; %take the RGB values
    % 
    % MapAnsColor = MapAnsColorTable(1:2:359, :); 
    %     %% Matlab can not use color map more than 255 value (?)
    %     %% so, use half (180) of the color table
    % figure; colormap(MapAnsColor);   
    % rgbplot(MapAnsColor); axis([0 180 -1 2]);
    % colorbar; 
    %     %% for some reaon the color table is not centered at blue
    %     %% so, I shifted it 
    % MapAnsColor= [MapAnsColor(31:180, :); MapAnsColor(1:30, :)];
    % figure; colormap(MapAnsColor); 
    % rgbplot(MapAnsColor); axis([0 180 -1 2]);
    % colorbar;
    %%%% saved in E:\Retinotopy\MapAnsColorWorkspace


% --- Outputs from this function are returned to the command line.
function varargout = ImanMapAns_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function TemplateValue_Callback(hObject, eventdata, handles)
% hObject    handle to TemplateValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TemplateValue as text
%        str2double(get(hObject,'String')) returns contents of TemplateValue as a double

global TemplateMethod TemplateValue;
TemplateValue = str2num(get(hObject, 'string'));

% --- Executes during object creation, after setting all properties.
function TemplateValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TemplateValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filter_size_edit_Callback(hObject, eventdata, handles)
% hObject    handle to filter_size_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filter_size_edit as text
%        str2double(get(hObject,'String')) returns contents of filter_size_edit as a double
global filter_size;
filter_size = str2num(get(hObject, 'string'));


% --- Executes during object creation, after setting all properties.
function filter_size_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_size_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LinePlotPushBtn.
function LinePlotPushBtn_Callback(hObject, eventdata, handles)
% hObject    handle to LinePlotPushBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

IMAPhase = IMA_draw('line');

% --- Executes on button press in SavePushBtn.
function SavePushBtn_Callback(hObject, eventdata, handles)
% hObject    handle to SavePushBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in RespAreaPushBtn.
function RespAreaPushBtn_Callback(hObject, eventdata, handles)
% hObject    handle to RespAreaPushBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global TemplateMethod TemplateValue;
if (TemplateMethod~=2&TemplateMethod~=3) 
    errordlg('wrong method for selecting threshold');
else
    IMA_measure;
end

% --- Executes on button press in ComparePhasePushBtn.
function ComparePhasePushBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ComparePhasePushBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
map2comp('phase_comp');
   

% --- Executes on button press in QuitPushBtn.
function QuitPushBtn_Callback(hObject, eventdata, handles)
% hObject    handle to QuitPushBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 clear all;
 close all; 

% --- Executes on button press in CompareAmpPushBtn.
function CompareAmpPushBtn_Callback(hObject, eventdata, handles)
% hObject    handle to CompareAmpPushBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
map2comp('amp_comp');


function spatial_freq_edit_Callback(hObject, eventdata, handles)
% hObject    handle to spatial_freq_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spatial_freq_edit as text
%        str2double(get(hObject,'String')) returns contents of spatial_freq_edit as a double
global deg_per_stimcycle;
deg_per_stimcycle = str2num(get(hObject, 'string'));

% --- Executes during object creation, after setting all properties.
function spatial_freq_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spatial_freq_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function um_per_pixel_edit_Callback(hObject, eventdata, handles)
% hObject    handle to um_per_pixel_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of um_per_pixel_edit as text
%        str2double(get(hObject,'String')) returns contents of um_per_pixel_edit as a double
global mm_per_pixel;
mm_per_pixel = str2num(get(hObject, 'string'))/1000;

% --- Executes during object creation, after setting all properties.
function um_per_pixel_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to um_per_pixel_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ContourPlotPushBtn.
function ContourPlotPushBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ContourPlotPushBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

IMAPhase = IMA_draw('contour');


function contourvalueeditnew_Callback(hObject, eventdata, handles)
% hObject    handle to contourvalueeditnew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of contourvalueeditnew as text
%        str2double(get(hObject,'String')) returns contents of contourvalueeditnew as a double


% --- Executes during object creation, after setting all properties.
function contourvalueeditnew_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contourvalueeditnew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MapQPushBtn.
function MapQPushBtn_Callback(hObject, eventdata, handles)
% hObject    handle to MapQPushBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
mapQ;

% --- Executes on selection change in TemplateMethod.
function TemplateMethod_Callback(hObject, eventdata, handles)
% hObject    handle to TemplateMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns TemplateMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TemplateMethod
global TemplateMethod TemplateValue;

TemplateMethod = get(hObject,'Value');
switch TemplateMethod
    case 1 %% exsited template
        set(handles.TemplateValue,'Enable','off');
    case 2 %% X fraction of the peak response
        set(handles.TemplateValue,'Enable','on');
        set(handles.TemplateValue,'string', '0.3'); %default value
    case 3 %% X fold of background
        set(handles.TemplateValue,'Enable','on');
        set(handles.TemplateValue,'string', '1.5'); %default value
    case 4 %% X pixels with strongest response
        set(handles.TemplateValue,'Enable','on');
        set(handles.TemplateValue,'string', '20000'); %default value
    case 5 %%read in a tif file generated by IDL
        set(handles.TemplateValue,'Enable','off');
    case 6 %% input local cordinates
        set(handles.TemplateValue,'Enable','off');
end
     
TemplateValue = str2num(get(handles.TemplateValue, 'string'));



% --- Executes during object creation, after setting all properties.
function TemplateMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TemplateMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in HelpPushBtn.
function HelpPushBtn_Callback(hObject, eventdata, handles)
% hObject    handle to HelpPushBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

IMA_help;
