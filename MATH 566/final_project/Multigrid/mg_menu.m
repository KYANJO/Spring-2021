function varargout = mg_menu(varargin)
% MG_MENU M-file for mg_menu.fig
%      MG_MENU, by itself, creates a new MG_MENU or raises the existing
%      singleton*.
%
%      H = MG_MENU returns the handle to a new MG_MENU or the handle to
%      the existing singleton*.
%
%      MG_MENU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MG_MENU.M with the given input arguments.
%
%      MG_MENU('Property','Value',...) creates a new MG_MENU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mg_menu_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mg_menu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mg_menu

% Last Modified by GUIDE v2.5 29-Mar-2004 18:30:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mg_menu_OpeningFcn, ...
                   'gui_OutputFcn',  @mg_menu_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mg_menu is made visible.
function mg_menu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mg_menu (see VARARGIN)

% Choose default command line output for mg_menu
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mg_menu wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mg_menu_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function g_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to g_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor','white');
end



function g_edit_Callback(hObject, eventdata, handles)
% hObject    handle to g_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of g_edit as text
%        str2double(get(hObject,'String')) returns contents of g_edit as a double


% --- Executes during object creation, after setting all properties.
function tol_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tol_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor','white');
end



function tol_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tol_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tol_edit as text
%        str2double(get(hObject,'String')) returns contents of tol_edit as a double


% --- Executes during object creation, after setting all properties.
function maxit_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxit_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor','white');
end



function maxit_edit_Callback(hObject, eventdata, handles)
% hObject    handle to maxit_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxit_edit as text
%        str2double(get(hObject,'String')) returns contents of maxit_edit as a double


% --- Executes during object creation, after setting all properties.
function cycle_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cycle_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cycle_popup.
function cycle_popup_Callback(hObject, eventdata, handles)
% hObject    handle to cycle_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns cycle_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cycle_popup


% --- Executes during object creation, after setting all properties.
function smoother_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smoother_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in smoother_popup.
function smoother_popup_Callback(hObject, eventdata, handles)
% hObject    handle to smoother_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns smoother_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from smoother_popup


% --- Executes during object creation, after setting all properties.
function omega_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to omega_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor','white');
end



function omega_edit_Callback(hObject, eventdata, handles)
% hObject    handle to omega_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of omega_edit as text
%        str2double(get(hObject,'String')) returns contents of omega_edit as a double


% --- Executes during object creation, after setting all properties.
function nu1_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nu1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor','white');
end



function nu1_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nu1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nu1_edit as text
%        str2double(get(hObject,'String')) returns contents of nu1_edit as a double


% --- Executes during object creation, after setting all properties.
function nu2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nu2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor','white');
end



function nu2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nu2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nu2_edit as text
%        str2double(get(hObject,'String')) returns contents of nu2_edit as a double


% --- Executes during object creation, after setting all properties.
function steps_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to steps_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor','white');
end



function steps_edit_Callback(hObject, eventdata, handles)
% hObject    handle to steps_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of steps_edit as text
%        str2double(get(hObject,'String')) returns contents of steps_edit as a double


% --- Executes during object creation, after setting all properties.
function Hmax_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Hmax_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor','white');
end



function Hmax_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Hmax_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Hmax_edit as text
%        str2double(get(hObject,'String')) returns contents of Hmax_edit as a double


% --- Executes on button press in exact_check.
function exact_check_Callback(hObject, eventdata, handles)
% hObject    handle to exact_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of exact_check


% --- Executes during object creation, after setting all properties.
function u_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to u_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor','white');
end



function u_edit_Callback(hObject, eventdata, handles)
% hObject    handle to u_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of u_edit as text
%        str2double(get(hObject,'String')) returns contents of u_edit as a double


% --- Executes during object creation, after setting all properties.
function b_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to b_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor','white');
end



function b_edit_Callback(hObject, eventdata, handles)
% hObject    handle to b_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of b_edit as text
%        str2double(get(hObject,'String')) returns contents of b_edit as a double


% --- Executes during object creation, after setting all properties.
function a_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor','white');
end



function a_edit_Callback(hObject, eventdata, handles)
% hObject    handle to a_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a_edit as text
%        str2double(get(hObject,'String')) returns contents of a_edit as a double


% --- Executes during object creation, after setting all properties.
function c_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor','white');
end



function c_edit_Callback(hObject, eventdata, handles)
% hObject    handle to c_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c_edit as text
%        str2double(get(hObject,'String')) returns contents of c_edit as a double


% --- Executes during object creation, after setting all properties.
function f_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor','white');
end



function f_edit_Callback(hObject, eventdata, handles)
% hObject    handle to f_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f_edit as text
%        str2double(get(hObject,'String')) returns contents of f_edit as a double


% --- Executes on button press in prec_check.
function prec_check_Callback(hObject, eventdata, handles)
% hObject    handle to prec_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of prec_check


% --- Executes on button press in direct_check.
function direct_check_Callback(hObject, eventdata, handles)
% hObject    handle to direct_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of direct_check


% --- Executes on button press in solve_push.
function solve_push_Callback(hObject, eventdata, handles)
% hObject    handle to solve_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
steps = str2double(get(handles.steps_edit,'String'));
Hmax = str2double(get(handles.Hmax_edit,'String'));
all_choices = get(handles.smoother_popup, 'String');
smooth = all_choices{get(handles.smoother_popup, 'Value')};
omega = str2double(get(handles.omega_edit,'String'));
nu1 = str2double(get(handles.nu1_edit,'String'));
nu2 = str2double(get(handles.nu2_edit,'String'));
all_choices = get(handles.cycle_popup, 'String');
cycle = all_choices{get(handles.cycle_popup, 'Value')};
if strcmp(cycle,'V-cycle') mu = 1; else mu = 2; end
tol = str2double(get(handles.tol_edit,'String'));
maxit = str2double(get(handles.maxit_edit,'String'));
geom = get(handles.g_edit,'String');
b = get(handles.b_edit,'String');
a = str2num(get(handles.a_edit,'String'));
c = str2double(get(handles.c_edit,'String'));
f = get(handles.f_edit,'String');
u = get(handles.u_edit,'String');
if (get(handles.prec_check, 'Value') == 1) prec = 1; else prec = 0; end
if (get(handles.exact_check, 'Value') == 1) exact = 1; else exact = 0; end   
if (get(handles.direct_check, 'Value') == 1) direct = 1; else direct = 0; end   
[p, e, t, u_lmax] = mg_main(steps, Hmax, smooth, omega, nu1, nu2, mu, tol, maxit, ...
			   exact, direct, prec, geom, b, a, c, f, u);
if (get(handles.plot_check, 'Value') == 1) 
  figure(1);
  pdesurf(p, t, u_lmax);
end
return;


% --- Executes on button press in plot_check.
function plot_check_Callback(hObject, eventdata, handles)
% hObject    handle to plot_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_check


