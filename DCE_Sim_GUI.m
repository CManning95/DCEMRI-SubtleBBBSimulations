function varargout = DCE_Sim_GUI(varargin)
% DCE_SIM_GUI MATLAB code for DCE_Sim_GUI.fig

% Will assess various physiological and technical DCE parameters for PS
% and vP accuracy

% Last Modified by GUIDE v2.5 21-Jan-2020 16:33:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DCE_Sim_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @DCE_Sim_GUI_OutputFcn, ...
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

addpath('DCE_Simulation_Functions');

% --- Executes just before DCE_Sim_GUI is made visible.
function DCE_Sim_GUI_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for DCE_Sim_GUI
handles.output = hObject;

% set default for plot_hold to 0
handles.plot_hold = 0;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = DCE_Sim_GUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function Hct_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Hct_Callback(hObject, eventdata, handles)
Hct = str2double(get(hObject,'String'));
if isnan(Hct)
    set(hObject, 'string', 0.45)
end
handles.PhysParam.Hct = Hct;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function vE_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function vE_Callback(hObject, eventdata, handles)
vE = str2double(get(hObject,'String'));
if isnan(vE)
    set(hObject, 'string', 0.2)
end
handles.PhysParam.vE = vE;
guidata(hObject, handles);

% ------------------------------------------------------------------------
function initialise_DCE_NAWM_MSS3(handles)

handles.PhysParam.Hct = 0.42;
handles.PhysParam.vE = 0.2;
handles.PhysParam.FP_mlPer100gPerMin = 18;
handles.PhysParam.T10_blood_s = 1.901;
handles.PhysParam.T10_tissue_s = 0.917;
handles.PhysParam.T2s0_blood_s = 0.191;
handles.PhysParam.T2s0_tissue_s = 0.050;
handles.PhysParam.S0_blood = 100;
handles.PhysParam.S0_tissue = 100;
handles.PhysParam.kbe_perS = 1.57;
handles.PhysParam.kie_perS = 1.81;
handles.PhysParam.vP_fixed = 0.0045;
handles.PhysParam.vP_fixed_single = 0.0045;
handles.PhysParam.PS_fixed = 1.5 * 1e-4;
handles.PhysParam.PS_fixed_single = 1.5 * 1e-4;

set(handles.Hct, 'String', handles.PhysParam.Hct);
set(handles.vE, 'String', handles.PhysParam.vE);
set(handles.Fp, 'String', handles.PhysParam.FP_mlPer100gPerMin);
set(handles.T10_blood_s, 'String', handles.PhysParam.T10_blood_s);
set(handles.T10_tissue_s, 'String', handles.PhysParam.T10_tissue_s);
set(handles.T2s0_blood_s, 'String', handles.PhysParam.T2s0_blood_s);
set(handles.T2s0_tissue_s, 'String',handles.PhysParam.T2s0_tissue_s);
set(handles.S0_blood, 'String', handles.PhysParam.S0_blood);
set(handles.S0_tissue, 'String', handles.PhysParam.S0_tissue);
set(handles.Kbe, 'String', handles.PhysParam.kbe_perS);
set(handles.Kie, 'String', handles.PhysParam.kie_perS);
set(handles.vP_fixed,'String',handles.PhysParam.vP_fixed);
set(handles.vP_fixed_single,'String',handles.PhysParam.vP_fixed_single);
set(handles.PS_fixed,'String',handles.PhysParam.PS_fixed * 1e4);
set(handles.PS_fixed_single,'String',handles.PhysParam.PS_fixed_single * 1e4);

handles.SeqParam.t_acq_s = 1268;
handles.SeqParam.t_res_sample_s = 39.62;
handles.SeqParam.TR_s = 1e-3*3.4;
handles.SeqParam.TE_s = 1e-3*1.7;
handles.SeqParam.r1_per_mM_per_s = 5.0;
handles.SeqParam.r2_per_mM_per_s = 7.1;
handles.SeqParam.FA_nom_deg = 15;
handles.SeqParam.FA_error = 1;

set(handles.t_acq_s, 'String', handles.SeqParam.t_acq_s);
set(handles.t_res_sample, 'String', handles.SeqParam.t_res_sample_s);
set(handles.TR,'String',1e3*handles.SeqParam.TR_s);
set(handles.TE,'String',1e3*handles.SeqParam.TE_s);
set(handles.r1,'String',handles.SeqParam.r1_per_mM_per_s);
set(handles.r2,'String',handles.SeqParam.r2_per_mM_per_s);
set(handles.FA_nom_deg,'String',handles.SeqParam.FA_nom_deg);
set(handles.FA_error,'String',handles.SeqParam.FA_error);

handles.SimParam.N_repetitions = 1000;
handles.SimParam.t_res_full_s = 0.1;
handles.SimParam.NIgnore = 0;
handles.SimParam.SNR = 164;
handles.SimParam.drift_pctPerMin = 0;
handles.SimParam.min_PS = 0 * 1e-4;
handles.SimParam.max_PS = 2 * 1e-4;
handles.SimParam.min_vP = 0;
handles.SimParam.max_vP = 0.01;
handles.SimParam.venous_delay_s = 6;
t_start = 0;
handles.SimParam.t_start_s = t_start*handles.SeqParam.t_res_sample_s;

set(handles.N_repetitions,'string',handles.SimParam.N_repetitions);
set(handles.t_res_full_s, 'string',handles.SimParam.t_res_full_s);
set(handles.NIgnore, 'string',handles.SimParam.NIgnore);
set(handles.SNR, 'string', handles.SimParam.SNR);
set(handles.drift_pctPerMin,'string', handles.SimParam.drift_pctPerMin);
set(handles.min_PS,'string', handles.SimParam.min_PS * 1e4);
set(handles.max_PS,'string', handles.SimParam.max_PS * 1e4);
set(handles.min_vP,'String', handles.SimParam.min_vP);
set(handles.max_vP,'String', handles.SimParam.max_vP);
set(handles.venous_delay_s,'string',handles.SimParam.venous_delay_s);
set(handles.t_start,'string',t_start);

handles.acqParam.T1_acq_method = 'Accurate';
handles.SimParam.InjectionRate = 'slow';
handles.SimParam.syn_model = '2CXM';
handles.SimParam.water_exch_model = 'FXL';
handles.acqParam.B1_correction = 'none';

set(handles.T1_acq_method,'value',3);
set(handles.InjectionRate,'value',2);
set(handles.syn_model, 'value',2);
set(handles.water_exch_model, 'value',2);
set(handles.B1_correction, 'value',0);

guidata(handles.figure1,handles)

function initialise_DCE_scGM_MSS3(handles)

handles.PhysParam.Hct = 0.42;
handles.PhysParam.vE = 0.2;
handles.PhysParam.FP_mlPer100gPerMin = 42;
handles.PhysParam.T10_blood_s = 1.901;
handles.PhysParam.T10_tissue_s = 1.212;
handles.PhysParam.T2s0_blood_s = 0.191;
handles.PhysParam.T2s0_tissue_s = 0.11;
handles.PhysParam.S0_blood = 100;
handles.PhysParam.S0_tissue = 100;
handles.PhysParam.kbe_perS = 1.83;
handles.PhysParam.kie_perS = 1.81;
handles.PhysParam.vP_fixed = 0.0127;
handles.PhysParam.vP_fixed_single = 0.0127;
handles.PhysParam.PS_fixed = 2.5 * 1e-4;
handles.PhysParam.PS_fixed_single = 2.5 * 1e-4;

set(handles.Hct, 'String', handles.PhysParam.Hct);
set(handles.vE, 'String', handles.PhysParam.vE);
set(handles.Fp, 'String', handles.PhysParam.FP_mlPer100gPerMin);
set(handles.T10_blood_s, 'String', handles.PhysParam.T10_blood_s);
set(handles.T10_tissue_s, 'String', handles.PhysParam.T10_tissue_s);
set(handles.T2s0_blood_s, 'String', handles.PhysParam.T2s0_blood_s);
set(handles.T2s0_tissue_s, 'String',handles.PhysParam.T2s0_tissue_s);
set(handles.S0_blood, 'String', handles.PhysParam.S0_blood);
set(handles.S0_tissue, 'String', handles.PhysParam.S0_tissue);
set(handles.Kbe, 'String', handles.PhysParam.kbe_perS);
set(handles.Kie, 'String', handles.PhysParam.kie_perS);
set(handles.vP_fixed,'String',handles.PhysParam.vP_fixed);
set(handles.vP_fixed_single,'String',handles.PhysParam.vP_fixed_single);
set(handles.PS_fixed,'String',handles.PhysParam.PS_fixed * 1e4);
set(handles.PS_fixed_single,'String',handles.PhysParam.PS_fixed_single * 1e4);

handles.SeqParam.t_acq_s = 1268;
handles.SeqParam.t_res_sample_s = 39.62;
handles.SeqParam.TR_s = 1e-3*3.4;
handles.SeqParam.TE_s = 1e-3*1.7;
handles.SeqParam.r1_per_mM_per_s = 5.0;
handles.SeqParam.r2_per_mM_per_s = 7.1;
handles.SeqParam.FA_nom_deg = 15;
handles.SeqParam.FA_error = 1;

set(handles.t_acq_s, 'String', handles.SeqParam.t_acq_s);
set(handles.t_res_sample, 'String', handles.SeqParam.t_res_sample_s);
set(handles.TR,'String',1e3*handles.SeqParam.TR_s);
set(handles.TE,'String',1e3*handles.SeqParam.TE_s);
set(handles.r1,'String',handles.SeqParam.r1_per_mM_per_s);
set(handles.r2,'String',handles.SeqParam.r2_per_mM_per_s);
set(handles.FA_nom_deg,'String',handles.SeqParam.FA_nom_deg);
set(handles.FA_error,'String',handles.SeqParam.FA_error);

handles.SimParam.N_repetitions = 1000;
handles.SimParam.t_res_full_s = 0.1;
handles.SimParam.NIgnore = 0;
handles.SimParam.SNR = 90;
handles.SimParam.drift_pctPerMin = 0;
handles.SimParam.min_PS = 0 * 1e-4;
handles.SimParam.max_PS = 5 * 1e-4;
handles.SimParam.min_vP = 0;
handles.SimParam.max_vP = 0.02;
handles.SimParam.venous_delay_s = 6;
t_start = 0;
handles.SimParam.t_start_s = t_start*handles.SeqParam.t_res_sample_s;

set(handles.N_repetitions,'string',handles.SimParam.N_repetitions);
set(handles.t_res_full_s, 'string',handles.SimParam.t_res_full_s);
set(handles.NIgnore, 'string',handles.SimParam.NIgnore);
set(handles.SNR, 'string', handles.SimParam.SNR);
set(handles.drift_pctPerMin,'string', handles.SimParam.drift_pctPerMin);
set(handles.min_PS,'string', handles.SimParam.min_PS * 1e4);
set(handles.max_PS,'string', handles.SimParam.max_PS * 1e4);
set(handles.min_vP,'String', handles.SimParam.min_vP);
set(handles.max_vP,'String', handles.SimParam.max_vP);
set(handles.venous_delay_s,'string',handles.SimParam.venous_delay_s);
set(handles.t_start,'string',t_start);

handles.acqParam.T1_acq_method = 'Accurate';
handles.SimParam.InjectionRate = 'slow';
handles.SimParam.syn_model = '2CXM';
handles.SimParam.water_exch_model = 'FXL';
handles.acqParam.B1_correction = 'none';

set(handles.T1_acq_method,'value',3);
set(handles.InjectionRate,'value',2);
set(handles.syn_model, 'value',2);
set(handles.water_exch_model, 'value',2);
set(handles.B1_correction, 'value',0);

guidata(handles.figure1,handles)

% --- Executes on button press in pushbutton1 - resets values to default
function pushbutton1_Callback(hObject, eventdata, handles)
initialise_DCE_NAWM_MSS3(handles);

function Fp_Callback(hObject, eventdata, handles)
Fp = str2double(get(hObject,'String'));
if isnan(Fp)
    set(hObject, 'string', 10);
end
handles.PhysParam.FP_mlPer100gPerMin = Fp;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Fp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function T10_blood_s_Callback(hObject, eventdata, handles)
T10_blood_s = str2double(get(hObject,'String'));
if isnan(T10_blood_s)
    set(hObject, 'string', 1.901);
end
handles.PhysParam.T10_blood = T10_blood_s;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function T10_blood_s_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function T10_tissue_s_Callback(hObject, eventdata, handles)
T10_tissue_s = str2double(get(hObject,'String'));
if isnan(T10_tissue_s);
    set(hObject, 'string', 0.917);
end
handles.PhysParam.T10_tissue_s = T10_tissue_s;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function T10_tissue_s_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function T2s0_blood_s_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function T2s0_blood_s_Callback(hObject, eventdata, handles)
T2s0_blood_s = str2double(get(hObject,'String'));
if isnan(T2s0_blood_s);
    set(hObject,'string',0.191);
end
handles.PhysParam.T2s0_blood = T2s0_blood_s;
guidata(hObject, handles);

function S0_blood_Callback(hObject, eventdata, handles)
S0_blood = str2double(get(hObject,'String'));
if isnan(S0_blood);
    set(hObject,'string',100);
end
handles.PhysParam.S0_blood = S0_blood;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function S0_blood_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function S0_tissue_Callback(hObject, eventdata, handles)
S0_tissue = str2double(get(hObject,'String'));
if isnan(S0_tissue);
    set(hObject,'string',100);
end
handles.PhysParam.S0_tissue = S0_tissue;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function S0_tissue_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Kbe_Callback(hObject, eventdata, handles)
Kbe = str2double(get(hObject,'String'));
if isnan(Kbe);
    set(hObject, 'string', 1.1);
end
handles.PhysParam.kbe_perS = Kbe;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Kbe_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Kie_Callback(hObject, eventdata, handles)
Kie = str2double(get(hObject,'String'));
if isnan(Kie);
    set(hObject,'string',1.81);
end
handles.PhysParam.kie_perS = Kie;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Kie_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in T1_acq_method.
function T1_acq_method_Callback(hObject, eventdata, handles)
T1_opt = get(handles.T1_acq_method, 'Value');
switch T1_opt
    case 1
    case 2
        handles.acqParam.T1_acq_method = 'VFA';
    case 3
        handles.acqParam.T1_acq_method = 'Accurate';
    case 4
        handles.acqParam.T1_acq_method = 'Assumed';
end
guidata(handles.figure1, handles);

function t_acq_s_Callback(hObject, eventdata, handles)
t_acq_s = str2double(get(hObject,'String'));
if isnan(t_acq_s)
    set(hObject,'string',1268);
end
handles.SeqParam.t_acq_s = t_acq_s;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function t_acq_s_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function t_res_sample_Callback(hObject, eventdata, handles)
t_res_sample_s = str2double(get(hObject,'String'));
if isnan(t_res_sample_s)
    set(hObject,'string',39.62);
end
handles.SeqParam.t_res_sample_s = t_res_sample_s;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function t_res_sample_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TR_Callback(hObject, eventdata, handles)
TR = str2double(get(hObject,'String'));
if isnan(TR)
    set(hObject,'string',3.4)
end
handles.SeqParam.TR_s = 1e-3 * TR;
guidata(hObject,handles)
    
% --- Executes during object creation, after setting all properties.
function TR_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TE_Callback(hObject, eventdata, handles)
TE = str2double(get(hObject,'String'));
if isnan(TE)
    set(hObject, 'string',1.7)
end
handles.SeqParam.TE_s = 1e-3 * TE;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function TE_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function r1_Callback(hObject, eventdata, handles)
r1 = str2double(get(hObject,'String'));
if isnan(r1)
    set(hObject,'string',5.0)
end
handles.SeqParam.r1_per_mM_per_s = r1;
guidata(hObject,handles)
    
% --- Executes during object creation, after setting all properties.
function r1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function r2_Callback(hObject, eventdata, handles)
 r2 = str2double(get(hObject,'String'));
 if isnan(r2)
     set(hObject,'string',7.1)
 end
 handles.SeqParam.r2_per_mM_per_s = r2;
 guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function r2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FA_nom_deg_Callback(hObject, eventdata, handles)
FA_nom_deg = str2double(get(hObject,'String'));
if isnan(FA_nom_deg)
    set(hObject,'string',15)
end
handles.SeqParam.FA_nom_deg = FA_nom_deg;
guidata(hObject,handles)
    
% --- Executes during object creation, after setting all properties.
function FA_nom_deg_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function N_repetitions_Callback(hObject, eventdata, handles)
N_repetitions = str2double(get(hObject,'String'));
if isnan(N_repetitions)
    set(hOject,'string',1000)
end
handles.SimParam.N_repetitions = N_repetitions;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function N_repetitions_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function t_res_full_s_Callback(hObject, eventdata, handles)
 t_res_full_s = str2double(get(hObject,'String'));
 if isnan(t_res_full_s)
     set(hObject, 'string', 0.1)
 end
 handles.SimParam.t_res_full_s = t_res_full_s;
 guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function t_res_full_s_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function NIgnore_Callback(hObject, eventdata, handles)
NIgnore = str2double(get(hObject,'String'));
if isnan(NIgnore) == 1 && strcmp(handles.SimParam.InjectionRate,'slow') == 1;
    set(hObject,'string',0)
elseif isnan(NIgnore) == 1 && strcmp(handles.SimParam.InjectionRate,'fast') == 1;
    set(hObject,'string',6)
end
    
handles.SimParam.NIgnore = 3 + NIgnore;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function NIgnore_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SNR_Callback(hObject, eventdata, handles)
SNR = str2double(get(hObject,'String'));
if isnan(SNR)
    set(hObject,'string',100)
end
handles.SimParam.SNR = SNR;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function SNR_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function drift_pctPerMin_Callback(hObject, eventdata, handles)
drift_pctPerMin = str2double(get(hObject,'String'));
if isnan(drift_pctPerMin)
    set(hObject,'string',0)
end
handles.SimParam.drift_pctPerMin = drift_pctPerMin;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function drift_pctPerMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FA_error_Callback(hObject, eventdata, handles)
 FA_error = str2double(get(hObject,'String'));
 if isnan(FA_error)
     set(hObject,'string',1)
 end
 handles.SeqParam.FA_error = FA_error;
 guidata(hObject,handles)
 
% --- Executes during object creation, after setting all properties.
function FA_error_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T2s0_tissue_s_Callback(hObject, eventdata, handles)
T2s0_tissue_s = str2double(get(hObject,'String'));
if isnan(T2s0_tissue_s)
    set(hObject,'string',0.080)
end
handles.PhysParam.T2s0_tissue_s = T2s0_tissue_s;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function T2s0_tissue_s_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Run_PS_Sims.
function Run_PS_Sims_Callback(hObject, eventdata, handles)
PhysParam = handles.PhysParam;
SeqParam = handles.SeqParam;
SimParam = handles.SimParam;

acqParam.T1_acq_method = handles.acqParam.T1_acq_method;
switch acqParam.T1_acq_method
    case 'VFA'
        acqParam.isFit = [0 0 1 1 1];
    case 'Accurate'
        acqParam.isFit = [1 1 1 1 1];
        acqParam.T1_acq_method = 'HIFI';
end

acqParam.TR_s = [0.0054 0.0054 0.0054 0.0054 0.0054];
acqParam.FA_nom_rads = [5 5 2 5 12] *2*(pi/360);
acqParam.FA_true_rads = SeqParam.FA_error * acqParam.FA_nom_rads;
acqParam.isIR = [1 1 0 0 0]; % indicates which are IR-SPGR
acqParam.TI_s = [0.168 1.068 NaN NaN NaN]; % Inversion times
acqParam.PECentre = [0.5 0.5 NaN NaN NaN]; % indicates time of centre of k-space
acqParam.NReadout = [160 160 160 160 160]; % number of readout pulses (Siemens - number of slices)
acqParam.NTry = 1;

if isempty('handles.assumed_T1_blood');
    handles.assumed_T1_blood = NaN;
end

if isempty('handles.assumed_T1_tissue');
    handles.assumed_T1_tissue = NaN;
end


% Simulate T1 acquisiton
[PhysParam.T1_blood_meas_s,temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,acqParam,acqParam.T1_acq_method,handles.assumed_T1_blood);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,acqParam,acqParam.T1_acq_method,handles.assumed_T1_tissue);

handles.FA_error_meas = acqParam.FA_error_meas;

disp(['Flip Angle error = ' num2str((100*SeqParam.FA_error)-100) ' %'])
disp(['Actual blood T1 = ' num2str(PhysParam.T10_blood_s)])
disp(['Measured blood T1 = ' num2str(PhysParam.T1_blood_meas_s)])
disp(['Actual tissue T1 = ' num2str(PhysParam.T10_tissue_s)])
disp(['Measured tissue T1 = ' num2str(PhysParam.T1_tissue_meas_s)])

switch handles.acqParam.B1_correction
    case 'HIFI'
        SeqParam.FA_meas_deg = acqParam.FA_error_meas * SeqParam.FA_nom_deg;
    case 'none'
        SeqParam.FA_meas_deg = SeqParam.FA_nom_deg;
end

% Check previous legend
if handles.plot_hold == 0; % delete previous figures if plot hold is off
    fig_GUI = findall(0,'type','figure','Tag','figure1'); %Keeps GUI open
    fig_other = findall(0,'type','figure');
    figures_del = setdiff(fig_other,fig_GUI);
    delete(figures_del);
end

if handles.plot_hold == 1 && ishandle(2) == 1; % if plot hold on, read current legend
    old_leg = (findall(figure(2),'Type','Legend'));
    old_legend.String = old_leg.String;
elseif handles.plot_hold == 1 && ishandle(2) == 0; % if no current legend, give actual PS string
    old_legend.String = ({'True PS'});
end

%Sort slow injection parameters
if SimParam.InjectionRate == 'slow'
    SimParam.baselineScans = 3; % datapoints to use for calculating base signal
    load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
    SimParam.Cp_AIF_mM = Cp_AIF_mM;
    SimParam.tRes_InputAIF_s = 18.49; % original time resolution of AIFs
    SimParam.InputAIFDCENFrames = 69; % number of time points
elseif SimParam.InjectionRate == 'fast'
    SimParam.baselineScans = 1; % datapoints to use for calculating base signal
end

%derive some additional parameters
SeqParam.NPoints = round(SeqParam.t_acq_s/SeqParam.t_res_sample_s);

SeqParam.FA_true_deg = SeqParam.FA_error*SeqParam.FA_meas_deg;

% Implement B1 correction if active
if handles.acqParam.B1_correction == 'HIFI'
    SeqParam.FA_meas_deg = SeqParam.FA_meas_deg*handles.FA_error_meas;
end

% ranges of PS to test for fixed vP
PS_range = linspace(SimParam.min_PS,SimParam.max_PS,10)+1e-8;
vP_fixed = PhysParam.vP_fixed; %vP for NAWM, Heye et al. 2016 

%range sizes
N_PS = size(PS_range,2);
PS_fit_1 = nan(SimParam.N_repetitions,N_PS);

for i_PS = 1:N_PS
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_1(:,i_PS)] = master_single_sim(PhysParam,SeqParam,SimParam);
end

PS_means_1 = mean(PS_fit_1,1); % mean for each PS
PS_devs_1 = std(PS_fit_1,0,1); % standard deviation for each PS

%% plotting graphs of results
figure(2);

set(gcf,'Units','centimeters','Position',[20,0,25,25]);

if exist('old_legend') == 1;
    if size(old_legend.String,2) == 1; % only plot true PS if first plot (only one entry in legend)
        plot( PS_range, zeros(size(PS_range)), 'k:','LineWidth',2); hold on;
    end
end

% errorbar(PS_range, PS_means_1(1,:), 1*PS_devs_1(1,:));
% %title(['Variable PS']);
% xlabel('PS (/min)');
% ylabel('fitted PS (/min)');
% xlim([0 SimParam.max_PS]);
% ylim([0*SimParam.max_PS SimParam.max_PS]);

errorbar(PS_range, PS_means_1(1,:) - PS_range, 1*PS_devs_1(1,:),'LineWidth',2);
%title(['Variable PS']);
xlabel('True PS (/min)');
ylabel('fitted PS error (/min)');
xlim([0 SimParam.max_PS]);
ylim([-3e-4 3e-4]);

if handles.plot_hold == 0; % if plot hold is off, overwrite legend
    legend({'True PS'}, 'Location', 'southeast');
    old_legend = findall(gcf,'Type','Legend');
end
if size(handles.plot_label,2) >= 1 % if theres a plot label, add it to legend
    old_legend.String{end+1} = handles.plot_label;
end

legend( old_legend.String , 'Location', 'northwest','FontSize',20 ); % plot updated legend

% set font size to 16 for axes ticks and labels
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',20);

% --- Executes on selection change in InjectionRate.
function InjectionRate_Callback(hObject, eventdata, handles)
Inj_Rate_opt = get(handles.InjectionRate,'Value');
switch Inj_Rate_opt
    case 1
    case 2
        handles.SimParam.InjectionRate = 'slow';
        handles.SimParam.NIgnore = 0;
        t_start = 0;
    case 3
        handles.SimParam.InjectionRate = 'fast';
        handles.SimParam.NIgnore = 6; % always ignore first 3 data points anyway (baseline)
        t_start = 198;
        %t_start = 3*handles.SeqParam.t_res_sample_s+5; % set min inj delay to 3 pre-contrast sample points(+5s)
end
set(handles.NIgnore,'string',handles.SimParam.NIgnore);
set(handles.t_start,'string',t_start);
handles.SimParam.t_start_s = t_start;
guidata(handles.figure1,handles)


% --- Executes during object creation, after setting all properties.
function InjectionRate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in syn_model.
function syn_model_Callback(hObject, eventdata, handles)
syn_model = get(handles.syn_model,'Value');
switch syn_model
    case 1
    case 2
        handles.SimParam.syn_model = '2CXM';
    case 3
        handles.SimParam.syn_model = 'Patlak';
end
guidata(handles.figure1,handles)



% --- Executes during object creation, after setting all properties.
function syn_model_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in water_exch_model.
function water_exch_model_Callback(hObject, eventdata, handles)
water_exch_model = get(handles.water_exch_model,'Value');
switch water_exch_model
    case 1
    case 2
        handles.SimParam.water_exch_model = 'FXL';
    case 3
        handles.SimParam.water_exch_model = 'SXL';
    case 4
        handles.SimParam.water_exch_model = '2S1XA';
    case 5
        handles.SimParam.water_exch_model = '3S2X_num';
end
guidata(handles.figure1,handles)

% --- Executes during object creation, after setting all properties.
function water_exch_model_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function min_PS_Callback(hObject, eventdata, handles)
min_PS = str2double(get(hObject,'String'));
if isnan(min_PS)
        set(hObject,'string',0)
end
handles.SimParam.min_PS = min_PS;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function min_PS_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function max_PS_Callback(hObject, eventdata, handles)
max_PS = str2double(get(hObject,'String'));
if isnan(max_PS)
    set(hObject,'string',5)
end
handles.SimParam.max_PS = max_PS * 1e-4;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function max_PS_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function vP_fixed_Callback(hObject, eventdata, handles)
vP_fixed = str2double(get(hObject,'String'));
if isnan(vP_fixed)
    set(hObject,'string',0.0058)
end
handles.PhysParam.vP_fixed = vP_fixed;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function vP_fixed_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function venous_delay_s_Callback(hObject, eventdata, handles)
venous_delay_s = str2double(get(hObject,'String'));
if isnan(venous_delay_s)
    set(hObject,'string',6)
end
handles.SimParam.venous_delay_s = venous_delay_s;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function venous_delay_s_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in B1_correction.
function B1_correction_Callback(hObject, eventdata, handles)
DCE_B1_correction = get(handles.B1_correction, 'Value');
switch DCE_B1_correction
    case 0
        handles.acqParam.B1_correction = 'none';
    case 1
        handles.acqParam.B1_correction = 'HIFI';
end
guidata(handles.figure1, handles);

% --- Executes during object creation, after setting all properties.
function B1_correction_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function t_start_Callback(hObject, eventdata, handles)
t_start = str2double(get(hObject,'String'));
if isnan(t_start) == 1 && strcmp(handles.SimParam.InjectionRate,'slow') == 1;
    set(hObject,'string',0)
elseif isnan(t_start) == 1 && strcmp(handles.SimParam.InjectionRate,'fast') == 1;
    set(hObject,'string',3*handles.SeqParam.t_res_sample_s + 5)
end

handles.SimParam.t_start_s = t_start;

if (strcmp(handles.SimParam.InjectionRate,'fast') == 1) && ((t_start < (3*handles.SeqParam.t_res_sample_s)) == 1)
    %error('Injection delay error')
    InjDelay_error = msgbox(['ERROR: Injection delay must be greater than 3 times the sample resolution, so that '...
        'there are at least 3 pre-contrast data points'], 'Error');
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function t_start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function min_vP_Callback(hObject, eventdata, handles)
min_vP = str2double(get(hObject,'String'));
if isnan(min_vP)
    set(hObject,'string',0)
end
handles.SimParam.min_vP = min_vP;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function min_vP_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_vP_Callback(hObject, eventdata, handles)
max_vP = str2double(get(hObject,'String'));
if isnan(max_vP)
    set(hObject,'string',0.02)
end
handles.SimParam.max_vP = max_vP;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function max_vP_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PS_fixed_Callback(hObject, eventdata, handles)
PS_fixed = str2double(get(hObject,'String'));
if isnan(PS_fixed)
    set(hObject,'string',2.5)
end
handles.PhysParam.PS_fixed = PS_fixed * 1e-4;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function PS_fixed_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run_vP_sims.
function run_vP_sims_Callback(hObject, eventdata, handles)
PhysParam = handles.PhysParam;
SeqParam = handles.SeqParam;
SimParam = handles.SimParam;

acqParam.T1_acq_method = handles.acqParam.T1_acq_method;
switch acqParam.T1_acq_method
    case 'VFA'
        acqParam.isFit = [0 0 1 1 1];
    case 'Accurate'
        acqParam.isFit = [1 1 1 1 1];
        acqParam.T1_acq_method = 'HIFI';
end

acqParam.TR_s = [0.0054 0.0054 0.0054 0.0054 0.0054];
acqParam.FA_nom_rads = [5 5 2 5 12] *2*(pi/360);
acqParam.FA_true_rads = SeqParam.FA_error * acqParam.FA_nom_rads;
acqParam.isIR = [1 1 0 0 0]; % indicates which are IR-SPGR
acqParam.TI_s = [0.168 1.068 NaN NaN NaN]; % Inversion times
acqParam.PECentre = [0.5 0.5 NaN NaN NaN]; % indicates time of centre of k-space
acqParam.NReadout = [160 160 160 160 160]; % number of readout pulses (Siemens - number of slices)
acqParam.NTry = 1;

if exist('handles.assumed_T1_blood') == 0;
    handles.assumed_T1_blood = NaN;
end

if exist('handles.assumed_T1_tissue') == 0;
    handles.assumed_T1_tissue = NaN;
end

% Simulate T1 acquisiton
[PhysParam.T1_blood_meas_s,temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,acqParam,acqParam.T1_acq_method,handles.assumed_T1_blood);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,acqParam,acqParam.T1_acq_method,handles.assumed_T1_tissue);

handles.FA_error_meas = acqParam.FA_error_meas;

disp(['Flip Angle error = ' num2str((100*SeqParam.FA_error)-100) ' %'])
disp(['Actual blood T1 = ' num2str(PhysParam.T10_blood_s)])
disp(['Measured blood T1 = ' num2str(PhysParam.T1_blood_meas_s)])
disp(['Actual tissue T1 = ' num2str(PhysParam.T10_tissue_s)])
disp(['Measured tissue T1 = ' num2str(PhysParam.T1_tissue_meas_s)])

switch handles.acqParam.B1_correction
    case 'HIFI'
        SeqParam.FA_meas_deg = acqParam.FA_error_meas * SeqParam.FA_nom_deg;
    case 'none'
        SeqParam.FA_meas_deg = SeqParam.FA_nom_deg;
end

% Check previous legend
if handles.plot_hold == 0; % delete previous figures if plot hold is off
    fig_GUI = findall(0,'type','figure','Tag','figure1'); %Keeps GUI open
    fig_other = findall(0,'type','figure');
    figures_del = setdiff(fig_other,fig_GUI);
    delete(figures_del);
end

if handles.plot_hold == 1 && ishandle(2) == 1; % if plot hold on, read current legend
    old_leg = (findall(figure(2),'Type','Legend'));
    old_legend.String = old_leg.String;
elseif handles.plot_hold == 1 && ishandle(2) == 0; % if no current legend, give actual PS string
    old_legend.String = ({'actual v_p'});
end

SimParam.baselineScans = 3; % datapoints to use for calculating base signal
%Sort slow injection parameters
if SimParam.InjectionRate == 'slow'
    load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
    SimParam.Cp_AIF_mM = Cp_AIF_mM;
    SimParam.tRes_InputAIF_s = 18.49; % original time resolution of AIFs
    SimParam.InputAIFDCENFrames = 69; % number of time points
end

%derive some additional parameters
SeqParam.NPoints = round(SeqParam.t_acq_s/SeqParam.t_res_sample_s);
SeqParam.FA_true_deg = SeqParam.FA_error*SeqParam.FA_meas_deg;

% Implement B1 correction if active
if handles.acqParam.B1_correction == 'HIFI'
    SeqParam.FA_meas_deg = SeqParam.FA_meas_deg*handles.FA_error_meas;
end

% ranges of PS to test for fixed vP
vP_range = linspace(SimParam.min_vP,SimParam.max_vP,10)+1e-8;
PS_fixed = PhysParam.PS_fixed; 

%range sizes
N_vP = size(vP_range,2);
vP_fit_1 = nan(SimParam.N_repetitions,N_vP);

for i_vP = 1:N_vP
    PhysParam.PS_perMin = PS_fixed(1);
    PhysParam.vP = vP_range(i_vP);
    [vP_fit_1(:,i_vP),temp] = master_single_sim(PhysParam,SeqParam,SimParam);
end

vP_means_1 = mean(vP_fit_1,1); % mean for each PS
vP_devs_1 = std(vP_fit_1,0,1); % standard deviation for each PS

%% plotting graphs of results
figure(2);
set(gcf,'Units','centimeters','Position',[20,0,20,20]);

% subplot(2,1,1)
if exist('old_legend') == 1;
    if size(old_legend.String,2) == 1; % only plot true vP if first plot (only one entry in legend)
        plot( vP_range, zeros(size(vP_range)), 'k:','LineWidth',2); hold on;
    end
end

% errorbar(vP_range, vP_means_1(1,:), 1*vP_devs_1(1,:)); hold on;
% title(['Variable vP']);
% xlabel('vP');
% ylabel('fitted vP');
% legend({'actual vP'}, 'Location', 'southeast');
% xlim([0 SimParam.max_vP]);
% ylim([0 SimParam.max_vP]);

errorbar(vP_range, vP_means_1(1,:) - vP_range, 1*vP_devs_1(1,:),'LineWidth',2);
xlabel('True vP');
ylabel('fitted vP error');
xlim([0 SimParam.max_vP]);
ylim([-7e-3 7e-3]);

if handles.plot_hold == 0; % if plot hold is off, overwrite legend
    legend({'actual vP'}, 'Location', 'southeast');
    old_legend = findall(gcf,'Type','Legend');
end
if size(handles.plot_label,2) >= 1 % if theres a plot label, add it to legend
    old_legend.String{end+1} = handles.plot_label;
end

legend( old_legend.String , 'Location', 'southeast' ); % plot updated legend
% set font size to 16 for axes ticks and labels
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',20);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
initialise_DCE_scGM_MSS3(handles);

% --- Executes on button press in plot_hold.
function plot_hold_Callback(hObject, eventdata, handles)
handles.plot_hold = get(hObject,'Value');
guidata(hObject,handles)

function plot_label_Callback(hObject, eventdata, handles)
handles.plot_label = get(hObject,'String');
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function plot_label_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text11_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in Single_PS.
function Single_PS_Callback(hObject, eventdata, handles)
PhysParam = handles.PhysParam;
SeqParam = handles.SeqParam;
SimParam = handles.SimParam;

acqParam.T1_acq_method = handles.acqParam.T1_acq_method;
switch acqParam.T1_acq_method
    case 'VFA'
        acqParam.isFit = [0 0 1 1 1];
    case 'Accurate'
        acqParam.isFit = [1 1 1 1 1];
        acqParam.T1_acq_method = 'HIFI';
end

acqParam.TR_s = [0.0054 0.0054 0.0054 0.0054 0.0054];
acqParam.FA_nom_rads = [5 5 2 5 12] *2*(pi/360);
acqParam.FA_true_rads = SeqParam.FA_error * acqParam.FA_nom_rads;
acqParam.isIR = [1 1 0 0 0]; % indicates which are IR-SPGR
acqParam.TI_s = [0.168 1.068 NaN NaN NaN]; % Inversion times
acqParam.PECentre = [0.5 0.5 NaN NaN NaN]; % indicates time of centre of k-space
acqParam.NReadout = [160 160 160 160 160]; % number of readout pulses (Siemens - number of slices)
acqParam.NTry = 1;

if exist('handles.assumed_T1_blood') == 0;
    handles.assumed_T1_blood = NaN;
end

if exist('handles.assumed_T1_tissue') == 0;
    handles.assumed_T1_tissue = NaN;
end


% Simulate T1 acquisiton
[PhysParam.T1_blood_meas_s,temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,acqParam,acqParam.T1_acq_method,handles.assumed_T1_blood);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,acqParam,acqParam.T1_acq_method,handles.assumed_T1_tissue);

handles.FA_error_meas = acqParam.FA_error_meas;

disp(['Flip Angle error = ' num2str((100*SeqParam.FA_error)-100) ' %'])
disp(['Actual blood T1 = ' num2str(PhysParam.T10_blood_s)])
disp(['Measured blood T1 = ' num2str(PhysParam.T1_blood_meas_s)])
disp(['Actual tissue T1 = ' num2str(PhysParam.T10_tissue_s)])
disp(['Measured tissue T1 = ' num2str(PhysParam.T1_tissue_meas_s)])

switch handles.acqParam.B1_correction
    case 'HIFI'
        SeqParam.FA_meas_deg = acqParam.FA_error_meas * SeqParam.FA_nom_deg;
    case 'none'
        SeqParam.FA_meas_deg = SeqParam.FA_nom_deg;
end


if handles.plot_hold == 0; % delete previous figures if plot hold is off
    fig_GUI = findall(0,'type','figure','Tag','figure1'); %Keeps GUI open
    fig_other = findall(0,'type','figure');
    figures_del = setdiff(fig_other,fig_GUI);
    delete(figures_del);
end

%Sort slow injection parameters
if SimParam.InjectionRate == 'slow'
    SimParam.baselineScans = 3; % datapoints to use for calculating base signal
    load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
    SimParam.Cp_AIF_mM = Cp_AIF_mM;
    SimParam.tRes_InputAIF_s = 18.49; % original time resolution of AIFs
    SimParam.InputAIFDCENFrames = 69; % number of time points
elseif SimParam.InjectionRate == 'fast'
    SimParam.baselineScans = 1; % datapoints to use for calculating base signal
end

%derive some additional parameters
SeqParam.NPoints = round(SeqParam.t_acq_s/SeqParam.t_res_sample_s);
SeqParam.FA_true_deg = SeqParam.FA_error*SeqParam.FA_meas_deg;

%Set PS value
PhysParam.vP = PhysParam.vP_fixed_single;
PhysParam.PS_perMin = PhysParam.PS_fixed_single;
[temp, PS_fit_1] = master_single_sim(PhysParam,SeqParam,SimParam);

CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'modal';
PS_msg = msgbox(['Mean PS = ' num2str(mean(PS_fit_1,2)*1e4) '(\pm ' num2str(1e4*1.96*std(PS_fit_1,0)) ') x 10^{-4} per min'], 'Single PS sim', CreateStruct);

function vP_fixed_single_Callback(hObject, eventdata, handles)
vP_fixed_single = str2double(get(hObject,'String'));
if isnan(vP_fixed_single)
    set(hObject,'string', 0.02)
end
handles.PhysParam.vP_fixed_single = vP_fixed_single;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function vP_fixed_single_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PS_fixed_single_Callback(hObject, eventdata, handles)
PS_fixed_single = str2double(get(hObject,'String'));
if isnan(PS_fixed_single)
    set(hObject,'string',2.5)
end
handles.PhysParam.PS_fixed_single = PS_fixed_single * 1e-4;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function PS_fixed_single_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function assumed_T1_blood_Callback(hObject, eventdata, handles)
assumed_T1_blood = str2double(get(hObject,'String'));
handles.assumed_T1_blood = assumed_T1_blood;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function assumed_T1_blood_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function assumed_T1_tissue_Callback(hObject, eventdata, handles)
assumed_T1_tissue = str2double(get(hObject,'String'));
handles.assumed_T1_tissue = assumed_T1_tissue;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function assumed_T1_tissue_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end