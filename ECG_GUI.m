function varargout = ECG_GUI(varargin)
% ECG_GUI MATLAB code for ECG_GUI.fig
%      ECG_GUI, by itself, creates a new ECG_GUI or raises the existing
%      singleton*.
%
%      H = ECG_GUI returns the handle to a new ECG_GUI or the handle to
%      the existing singleton*.
%
%      ECG_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ECG_GUI.M with the given input arguments.
%
%      ECG_GUI('Property','Value',...) creates a new ECG_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ECG_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ECG_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ECG_GUI

% Last Modified by GUIDE v2.5 07-Aug-2016 08:35:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ECG_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ECG_GUI_OutputFcn, ...
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


% --- Executes just before ECG_GUI is made visible.
function ECG_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ECG_GUI (see VARARGIN)

% Choose default command line output for ECG_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes ECG_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ECG_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load_ECG.
function load_ECG_Callback(hObject, eventdata, handles)
% hObject    handle to load_ECG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiopen('MAT')
grid on
xlabel('sample number')
ylabel('votage level')
a=val(1,:);
[vect,R_r,hrt,qrs_dur,beat,ECG]=analyzeECG(a);
plot(horzcat(ECG(1:360)))
set(handles.r_r,'String',R_r);
set(handles.heart_rate,'String',hrt);
set(handles.qrs_comp,'String',qrs_dur);
y=ECGNeuralNetworkFunction2(vect');
set(handles.result,'String',classify(y));



% --- Executes on button press in print.
function print_Callback(hObject, eventdata, handles)
% hObject    handle to print (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.result,'String','Printed');


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over print.
function [class_vect]=analyzeECG2(ecg_signal)   %for each beat
sig=ecg_signal;
n=ceil(log2(length(sig)));
f=360;
%%detrending
[p,s,mu] = polyfit((1:numel(sig)),sig,6);
f_y = polyval(p,(1:numel(sig)),[],mu);
ECG_data = sig-f_y;
%%analysis proper
f_comp=fft(ECG_data,2^n);
sig_energy=abs(f_comp.*f_comp);
%% my algorithm on peak finding
[locs_Qwave,locs_Rwave,~,qrs_area,qrs_duration,~]...
    =qrs_comp2(ECG_data,f);
%% heart rate calculation and QRS complex extraction
%END of QRS evaluation
%% T-wave detection
[~,t_peak,qt_interval,rt_interval,~,~]=t_wave2(ECG_data,f,locs_Qwave,locs_Rwave);
%% P-WAVE DETECTION
[~,pr_interval,p_peak,~,~,~]=p_wave2(ECG_data,f,locs_Qwave);
class_vect=[qrs_duration,qrs_area,pr_interval,rt_interval,qt_interval,t_peak,p_peak,sig_energy];


function heart_rate=HRT(locs_Rwave,f_sampling)
avg=[];
j=1;
for i=2:length(locs_Rwave)
    avg(j)=locs_Rwave(i)-locs_Rwave(i-1);
    j=j+1;
end
heart_rate=(sum(avg)/length(avg))*(1/f_sampling);

function [locs_Qwave,locs_Rwave,locs_Swave,qrs_area,qrs_duration,q_samples]...
    =qrs_comp(ECG_data,fs)
offset=round(0.2*fs);
%QRS COMPLEX
[R_amp,locs_Rwave] = findpeaks(ECG_data,'MinPeakHeight',max(ECG_data)*0.4,...
                                   'MinPeakDistance',100);
%for S wave
ECG_data=-ECG_data;
local_set=ECG_data(locs_Rwave:locs_Rwave+offset);
%local_set=-(local_set);
[~,ind]=findpeaks(local_set,'minpeakheight',0.4*max(local_set),'minpeakdistance',offset-1);
locs_Swave=locs_Rwave+ind-1;
if (locs_Swave(end)>length(ECG_data))
    locs_Swave(end)=length(ECG_data);
end
%end S wave
%for Q wave
if ((locs_Rwave(1)-offset)<=0)
    local_setQ=ECG_data(1:locs_Rwave);
else
    local_setQ=ECG_data((locs_Rwave-offset):locs_Rwave);
end
%local_setQ=-(local_setQ);
[~,indQ]=findpeaks(local_setQ,'minpeakheight',0.4*max(local_setQ),'minpeakdistance',offset-1);
locs_Qwave=locs_Rwave-(length(local_setQ)-indQ-1);
%%  qrs bla bla bla
qrs_duration=((locs_Swave-locs_Qwave)/fs)*1000;
qrs_duration=sum(qrs_duration)/length(qrs_duration);
qrs_area=0.5.*R_amp.*(qrs_duration/1000);
qrs_area=sum(qrs_area)/length(qrs_area);
%% for q_samples
sample_space=locs_Swave-locs_Qwave;
q_samples=ECG_data(locs_Swave:round(sample_space/6):locs_Qwave);
function [qt_int,t_sample,t_peak,qt_interval,rt_interval,sample_int,peak_intt]=t_wave(ECG_data,fs,locs_Qwave,locs_Rwave)
%  T-onset and off-set must be considered
%  the peak between duration=T_onset-T_offset is the T-wave
Q_T1=ceil(0.44*fs);
Q_T2=ceil(0.4*fs);
Q_Tm=ceil(0.36*fs);
% intervals
if (locs_Qwave(end)+Q_T1>fs*10)
    int1=locs_Qwave(1:end-1)+Q_T1;
    int2=locs_Qwave(1:end-1)+Q_T2;
    intm=locs_Qwave(1:end-1)+Q_Tm;
else
    int1=locs_Qwave(1:end)+Q_T1;
    int2=locs_Qwave(1:end)+Q_T2;
    intm=locs_Qwave(1:end)+Q_Tm;    
end
% T-wave interval
T_int=ceil(0.12*fs);
%sub_window
for i=1:length(int1)
    Q_section1(i,:)=ECG_data(int1(i)-T_int:int1(i));
    Q_section2(i,:)=ECG_data(int2(i)-T_int:int2(i));
    Q_sectionm(i,:)=ECG_data(intm(i)-T_int:intm(i));
end
%sub_peaks
[row,~]=size(Q_section1);
for i=1:row
    [T1(i),ind1(i)]=max(Q_section1(i,:));
    [T2(i),ind2(i)]=max(Q_section2(i,:));
    [Tm(i),indm(i)]=max(Q_sectionm(i,:));
end
t=[Tm; T1; T2];
% make  decision on peaks
[~,col]=size(t);
for i=1:col
    [t_peak(i),ind]=max(t(:,i));
    if ind==2
        qt_int(i,:)=(int1(i)-T_int):int1(i);
        qt_interval=(Q_T1/fs)*1000;
        RT_int(i)=int1(i)-locs_Rwave(i);
        rt_interval(i)=(RT_int(i)/fs)*1000;
        t_sample(i,:)=ECG_data(int1(i)-T_int:round(length(Q_section1)/6):int1(i));
        sample_int(i,:)=int1(i)-T_int:round(length(Q_section1)/6):int1(i);
        peak_intt(i)=(int1(i)-T_int)+ind1(i);
    elseif ind==3
        qt_int(i,:)=(int2(i)-T_int):int2(i);
        qt_interval=(Q_T2/fs)*1000;
        RT_int(i)=int2(i)-locs_Rwave(i);
        rt_interval(i)=(RT_int(i)/fs)*1000;
        t_sample(i,:)=ECG_data(int2(i)-T_int:round(length(Q_section2)/6):int2(i));
        sample_int(i,:)=int2(i)-T_int:round(length(Q_section2)/6):int2(i);
        peak_intt(i)=(int2(i)-T_int)+ind2(i);
    elseif ind==1
        qt_int(i,:)=(intm(i)-T_int):intm(i);
        qt_interval=(Q_Tm/fs)*1000;
        RT_int(i)=intm(i)-locs_Rwave(i);
        rt_interval(i)=(RT_int(i)/fs)*1000;
        t_sample(i,:)=ECG_data(intm(i)-T_int:round(length(Q_sectionm)/6):intm(i));
        sample_int(i,:)=intm(i)-T_int:round(length(Q_sectionm)/6):intm(i);
        peak_intt(i)=(intm(i)-T_int)+indm(i);
    end    
end
% t_peak=sum(t_peaks)/length(t_peaks);
rt_interval=sum(rt_interval)/length(rt_interval);

function [pr_int,pr_interval,p_peak,p_sample,p_int,peak_int]=p_wave(ECG_data,fs,locs_Qwave)
% declaring time range 5-6 small boxes
PQ_intSec1=4*0.04;
PQ_intSec2=5*0.04;
%range declared %convert them to samples
PQ_sample1=ceil(PQ_intSec1*fs);
PQ_sample2=ceil(PQ_intSec2*fs);
P_int=ceil(0.1*fs); %sample that comproses a P-wave
%sample conversion finish
origin_zero1=locs_Qwave-PQ_sample1;
origin_zero2=locs_Qwave-PQ_sample2;
if origin_zero1(1)<=0
    origin_zero1(1)=1;
    for i=2:length(locs_Qwave)
        P1_int(i,:)=origin_zero1(i):origin_zero1(i)+P_int;
    end
else
    for i=1:length(locs_Qwave)
        P1_int(i,:)=origin_zero1(i):origin_zero1(i)+P_int;
    end
end
if origin_zero2(1)<=0
    origin_zero2(1)=1;
    for i=2:length(locs_Qwave)
        P2_int(i,:)=origin_zero2(i):origin_zero2(i)+P_int;
    end
else
    for i=1:length(locs_Qwave)
        P2_int(i,:)=origin_zero2(i):origin_zero2(i)+P_int;
    end
end
[row,~]=size(P1_int);
P_section1=ECG_data(P1_int);
P_section2=ECG_data(P2_int);
% taking their peaks
for i=1:row
    [P1(i),ind1(i)]=max(P_section1(i,:));
    [P2(i),ind2(i)]=max(P_section2(i,:));
end
p=[P1;P2];
[~,col]=size(p);
%make decision peaks
for i=1:col
    [p_peak(i),ind]=max(p(:,i));
    if ind==1
        pr_int=P1_int;
        pr_interval=(length(pr_int(i,:))/fs)*1000;
        p_sample(i,:)=ECG_data(origin_zero1(i):round(length(pr_int(i,:))/6):pr_int(i,end));
        p_int(i,:)=origin_zero1(i):round(length(pr_int(i,:))/6):pr_int(i,end);
        peak_int(i)=ind1(i)+origin_zero1(i);
    elseif ind==2
        pr_int=P2_int;
        pr_interval=(length(pr_int(i,:))/fs)*1000;
        p_sample(i,:)=ECG_data(origin_zero2(i):round(length(pr_int(i,:))/6):pr_int(i,end));
        p_int(i,:)=origin_zero2(i):round(length(pr_int(i,:))/6):pr_int(i,end);
        peak_int(i)=ind2(i)+origin_zero2(i);
    end
end
% p_peak=sum(p_peak)/length(p_peak);
function resultString=classify(vect)
for i=1:4
    if vect(i)<0.5
        vect(i)=0;
    else
        vect(i)=1;
    end
end
if (vect(1) && not(vect(2)) && not(vect(3)) && not(vect(4)))
    resultString=sprintf('NORMAL');
elseif  (not(vect(1)) && vect(2) && not(vect(3)) && not(vect(4)))
    resultString=sprintf('LEFT BUNDLE BRANCH BLOCK');
elseif (not(vect(1)) && not(vect(2)) && vect(3) && not(vect(4)))
    resultString=sprintf('RIGHT BUNDLE BRANCH BLOCK');
elseif (not(vect(1)) && not(vect(2)) && not(vect(3)) && vect(4))
    resultString=sprintf('PACED');
else
    resultString=sprintf('OTHERS');
end
function [class_vect,R_r,hrt,qrs_dur,beat,ECG_data]=analyzeECG(ecg_signal) %for the whole record
sig=ecg_signal;
t=10;    %duration of signal
f=360;
%%detrending
[p,s,mu] = polyfit((1:numel(sig)),sig,6);
f_y = polyval(p,(1:numel(sig)),[],mu);
ECG_data = sig-f_y;
ECG_data=ECG_data(1:(f*t));
% end detrending
for i=1:10
    ECG_data=medfilt1(ECG_data,4);
end
% filtering end
%%analysis proper
%% my algorithm on peak finding
[locs_Qwave,locs_Rwave,~,~,qrs_duration,~]...
    =qrs_comp(ECG_data,f);
%% heart rate calculation and QRS complex extraction
heart_rate=HRT(locs_Rwave,f);
R_R=heart_rate*1000;
bpm=(1/heart_rate)*60;
R_r=sprintf(' %.0f ms',R_R);
hrt=sprintf(' %.0f bpm',bpm);
qrs_dur=sprintf(' %.0f ms',round(qrs_duration));
%END of QRS evaluation
%% T-wave detection
[~,~,~,~,~,~,locs_Twave]=t_wave(ECG_data,f,locs_Qwave,locs_Rwave);
%% P-WAVE DETECTION
[~,~,~,~,~,locs_Pwave]=p_wave(ECG_data,f,locs_Qwave);
%%classification vector
beat=beatExtractor(locs_Rwave,locs_Pwave,locs_Twave,ECG_data);
class_vect=analyzeECG2(beat(1,:));

function [qt_int,t_peak,qt_interval,rt_interval,sample_int,locs_Twave]=t_wave2(ECG_data,fs,locs_Qwave,locs_Rwave)
%  T-onset and off-set must be considered
%  the peak between duration=T_onset-T_offset is the T-wave
Q_T1=ceil(0.44*fs);
Q_T2=ceil(0.4*fs);
Q_Tm=ceil(0.36*fs);
% intervals
int1=locs_Qwave+Q_T1;
int2=locs_Qwave+Q_T2;
intm=locs_Qwave+Q_Tm;
% T-wave interval
T_int=ceil(0.12*fs);
%sub_window
if (int1>length(ECG_data))
    Q_section1(1,:)=ECG_data(int1-T_int:length(ECG_data));
    Q_section2(1,:)=ECG_data(int2-T_int:length(ECG_data));
    Q_sectionm(1,:)=ECG_data(intm-T_int:length(ECG_data));
else
    Q_section1(1,:)=ECG_data(int1-T_int:int1);
    Q_section2(1,:)=ECG_data(int2-T_int:int2);
    Q_sectionm(1,:)=ECG_data(intm-T_int:intm);
end
%sub_peaks
[T1,ind1]=max(Q_section1);
[T2,ind2]=max(Q_section2);
[Tm,indm]=max(Q_sectionm);
t=[Tm; T1; T2];
% make  decision on peaks
[t_peak,ind]=max(t);
if ind==2
  qt_int=(int1-T_int):int1;
  qt_interval=(Q_T1/fs)*1000;
  RT_int=int1-locs_Rwave;
  rt_interval=(RT_int/fs)*1000;
  %t_sample=ECG_data(int1-T_int:round(length(Q_section1)/6):int1);
  sample_int=int1-T_int:round(length(Q_section1)/6):int1;
  locs_Twave=(int1-T_int)+ind1;
elseif ind==3
  qt_int=(int2-T_int):int2;
  qt_interval=(Q_T2/fs)*1000;
  RT_int=int2-locs_Rwave;
  rt_interval=(RT_int/fs)*1000;
 % t_sample=ECG_data(int2-T_int:round(length(Q_section2)/6):int2);
  sample_int=int2-T_int:round(length(Q_section2)/6):int2;
  locs_Twave=(int2-T_int)+ind2;
elseif ind==1
  qt_int=(intm-T_int):intm;
  qt_interval=(Q_Tm/fs)*1000;
  RT_int=intm-locs_Rwave;
  rt_interval=(RT_int/fs)*1000;
 % t_sample=ECG_data(intm-T_int:round(length(Q_sectionm)/6):intm);
  sample_int=intm-T_int:round(length(Q_sectionm)/6):intm;
  locs_Twave=(intm-T_int)+indm;
end    
function [pr_int,pr_interval,p_peak,p_sample,p_int,peak_int]=p_wave2(ECG_data,fs,locs_Qwave)
% declaring time range 5-6 small boxes
PQ_intSec1=4*0.04;
PQ_intSec2=5*0.04;
%range declared %convert them to samples
PQ_sample1=ceil(PQ_intSec1*fs);
PQ_sample2=ceil(PQ_intSec2*fs);
P_int=ceil(0.1*fs); %sample that comproses a P-wave
%sample conversion finish
origin_zero1=locs_Qwave-PQ_sample1;
origin_zero2=locs_Qwave-PQ_sample2;
if origin_zero1<=0
    origin_zero1=1;
    P1_int=origin_zero1:origin_zero1+P_int;
else
    P1_int=origin_zero1:origin_zero1+P_int;
end
if origin_zero2<=0
    origin_zero2=1;
    P2_int=origin_zero2:origin_zero2+P_int;
else
    P2_int=origin_zero2:origin_zero2+P_int;
end
P_section1=ECG_data(P1_int);
P_section2=ECG_data(P2_int);
% taking their peaks
[P1,ind1]=max(P_section1);
ind1=origin_zero1+ind1;
[P2,ind2]=max(P_section2);
ind2=origin_zero2+ind2;
p=[P1;P2];
%make decision peaks
[p_peak,ind]=max(p);
if ind==1
     pr_int=P1_int;
     pr_interval=(length(pr_int)/fs)*1000;
     p_sample=ECG_data(origin_zero1:round(length(pr_int)/6):pr_int(1,end));
     p_int=origin_zero1:round(length(pr_int)/6):pr_int(1,end);
     peak_int=ind1;
elseif ind==2
     pr_int=P2_int;
     pr_interval=(length(pr_int)/fs)*1000;
     p_sample=ECG_data(origin_zero2:round(length(pr_int)/6):pr_int(1,end));
     p_int=origin_zero2:round(length(pr_int)/6):pr_int(1,end);
     peak_int=ind2;
end    
function beat=beatExtractor(locs_Rwave,locs_Pwave,locs_Twave,ECG_data)
%%evaluate positions
%locs_Pwave=locs_Pwave(:,1);
%locs_Twave=locs_Twave(:,end);
pos=[length(locs_Rwave),length(locs_Pwave),length(locs_Twave)];
iteration=min(pos);
max_length=max(locs_Twave(1:iteration)-locs_Pwave(1:iteration));
 for i=1:iteration
     step_length=locs_Twave(i)-locs_Pwave(i);
     beat(i,:)=horzcat(ECG_data(locs_Pwave(i):locs_Twave(i)),zeros(1,(max_length-step_length)));
 end
 function [locs_Qwave,locs_Rwave,locs_Swave,qrs_area,qrs_duration,q_samples]...
    =qrs_comp2(ECG_data,fs)
offset=round(0.2*fs);
%QRS COMPLEX
[R_amp,locs_Rwave] = findpeaks(ECG_data,'MinPeakHeight',max(ECG_data)*0.4,...
                                   'npeaks',1);
%for S wave
ECG_data=-ECG_data;
if (locs_Rwave+offset>length(ECG_data))
    local_set=ECG_data(locs_Rwave:end);
else
    local_set=ECG_data(locs_Rwave:locs_Rwave+offset);
end
%local_set=-(local_set);
[~,ind]=findpeaks(local_set,'minpeakheight',0.4*max(local_set),'npeaks',1);
locs_Swave=locs_Rwave+ind-1;
if (locs_Swave(end)>length(ECG_data))
    locs_Swave(end)=length(ECG_data);
end
%end S wave
%for Q wave
if ((locs_Rwave(1)-offset)<=0)
    local_setQ=ECG_data(1:locs_Rwave);
else
    local_setQ=ECG_data((locs_Rwave-offset):locs_Rwave);
end
%local_setQ=-(local_setQ);
[~,indQ]=findpeaks(local_setQ,'minpeakheight',0.4*max(local_setQ),'npeaks',1);
locs_Qwave=locs_Rwave-(length(local_setQ)-indQ-1);
%%  qrs bla bla bla
qrs_duration=((locs_Swave-locs_Qwave)/fs)*1000;
qrs_duration=sum(qrs_duration)/length(qrs_duration);
qrs_area=0.5.*R_amp.*(qrs_duration/1000);
qrs_area=sum(qrs_area)/length(qrs_area);
%% for q_samples
sample_space=locs_Swave-locs_Qwave;
q_samples=ECG_data(locs_Swave:round(sample_space/6):locs_Qwave);

function [Y,Xf,Af] = ECGNeuralNetworkFunction2(X,~,~)
%MYNEURALNETWORKFUNCTION neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 15-Aug-2016 06:05:02.
% 
% [Y] = myNeuralNetworkFunction(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timsteps
%   Each X{1,ts} = 263xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 4xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

  % ===== NEURAL NETWORK CONSTANTS =====
  
  % Input 1
  x1_step1_xoffset = [0;0;0;0;0;-2.48794001150818;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
  x1_step1_gain = [0.00791208791208791;0.0803540640625967;0.0194594594594595;0.00529411764705882;0.00452830188679245;0.0234541033165955;0.0374485269317441;2.5184106844152e-05;2.60401522828397e-05;1.04248200592073e-06;2.91063935155402e-08;2.31099261787283e-08;3.12060376399016e-08;6.01336096265736e-08;7.85323964703405e-08;9.81961469196636e-08;2.07371729705916e-07;2.35595084467373e-07;3.2455368184094e-07;4.85325949132636e-07;4.72183425894294e-07;6.90687553718566e-07;1.07944041418163e-06;1.65046916188026e-06;1.2977564605269e-06;2.28162995985367e-06;3.78133447092604e-06;3.79315283637669e-06;3.78981292567684e-06;7.82720984794236e-06;6.03300053032193e-06;1.25537167493668e-05;2.54438514834636e-05;2.47561451463498e-05;3.60689466839519e-05;5.21366475081503e-05;8.10060722297695e-05;0.000109211332466393;0.000149719481490066;0.000145214600514033;0.000159597792780454;0.000161441567253212;0.000165530148249129;0.000178897494998616;0.000193496607058083;0.00021245294159565;0.000227898498690591;0.000239350591074203;0.000255024738963992;0.00026917905275138;0.000272683454272654;0.000288474024791892;0.000300261918055845;0.000318208523686631;0.00034700700759093;0.00037327117594609;0.000395103530218371;0.00040961858460154;0.000417327130685508;0.000425651418007188;0.000442247980625904;0.00046275174249659;0.00049264464739158;0.000534746009120479;0.000569134865419573;0.000579490515822509;0.000580214183492575;0.000577305153014163;0.000587104261820675;0.000628639905586606;0.000683009365051442;0.000722253386048259;0.000745178885629688;0.000752128703573494;0.00074534872720981;0.000756210339001117;0.000783196197622161;0.00081672539297661;0.000861776762050294;0.000910220432547675;0.000929498026009277;0.000930009705222816;0.000930441639559274;0.000939782282705227;0.000975544329050371;0.00103191256319848;0.00107549195419326;0.00108968149864418;0.00109072583482341;0.00109014614237461;0.00110954829672921;0.00114648531541527;0.00119020835205863;0.00122940366031518;0.00125823271676889;0.00126627448344563;0.00126404858517551;0.00126071565075184;0.00127188872873045;0.0013173052772691;0.0013795687146583;0.0014256340019487;0.00144002978707101;0.00142609096898341;0.0014011547284368;0.00141040602995441;0.0014486175130365;0.00149079505098873;0.00153165979035303;0.00156385695024263;0.00156483805715054;0.00155864800431607;0.00155804581331634;0.00155503149709887;0.00157555355255662;0.00162593558216595;0.00166458717540643;0.00167183081878557;0.00166405010217779;0.00164124591997432;0.00163854813745586;0.00167399227598227;0.00171656057878072;0.00173562758631511;0.00174146952766769;0.00172722489148853;0.00170684299949076;0.00170370416918492;0.00171357982455891;0.00173314509547782;0.00176845764736203;0.00179905712146399;0.00179538949621084;0.00176454010033724;0.00172205719104508;0.00170283044047219;0.00172205719104508;0.00176454010033724;0.00179538949621084;0.00179905712146399;0.00176845764736203;0.00173314509547782;0.00171357982455891;0.00170370416918492;0.00170684299949076;0.00172722489148853;0.00174146952766769;0.00173562758631511;0.00171656057878072;0.00167399227598227;0.00163854813745586;0.00164124591997432;0.00166405010217779;0.00167183081878557;0.00166458717540643;0.00162593558216595;0.00157555355255662;0.00155503149709887;0.00155804581331634;0.00155864800431607;0.00156483805715054;0.00156385695024263;0.00153165979035303;0.00149079505098873;0.0014486175130365;0.00141040602995441;0.0014011547284368;0.00142609096898341;0.00144002978707101;0.0014256340019487;0.0013795687146583;0.0013173052772691;0.00127188872873045;0.00126071565075184;0.00126404858517551;0.00126627448344563;0.00125823271676889;0.00122940366031518;0.00119020835205863;0.00114648531541527;0.00110954829672921;0.00109014614237461;0.00109072583482341;0.00108968149864418;0.00107549195419326;0.00103191256319848;0.000975544329050371;0.000939782282705227;0.000930441639559274;0.000930009705222816;0.000929498026009277;0.000910220432547675;0.000861776762050294;0.00081672539297661;0.000783196197622161;0.000756210339001117;0.00074534872720981;0.000752128703573494;0.000745178885629688;0.000722253386048259;0.000683009365051442;0.000628639905586606;0.000587104261820675;0.000577305153014163;0.000580214183492575;0.000579490515822509;0.000569134865419573;0.000534746009120479;0.00049264464739158;0.00046275174249659;0.000442247980625904;0.000425651418007188;0.000417327130685508;0.00040961858460154;0.000395103530218371;0.00037327117594609;0.00034700700759093;0.000318208523686631;0.000300261918055845;0.000288474024791892;0.000272683454272654;0.00026917905275138;0.000255024738963992;0.000239350591074203;0.000227898498690591;0.00021245294159565;0.000193496607058083;0.000178897494998616;0.000165530148249129;0.000161441567253212;0.000159597792780454;0.000145214600514033;0.000149719481490066;0.000109211332466393;8.10060722297695e-05;5.21366475081503e-05;3.60689466839519e-05;2.47561451463498e-05;2.54438514834636e-05;1.25537167493668e-05;6.03300053032193e-06;7.82720984794236e-06;3.78981292567684e-06;3.79315283637669e-06;3.78133447092604e-06;2.28162995985367e-06;1.2977564605269e-06;1.65046916188026e-06;1.07944041418163e-06;6.90687553718566e-07;4.72183425894294e-07;4.85325949132636e-07;3.2455368184094e-07;2.35595084467373e-07;2.07371729705916e-07;9.81961469196636e-08;7.85323964703405e-08;6.01336096265736e-08;3.12060376399016e-08;2.31099261787283e-08;2.91063935155402e-08;1.04248200592073e-06;2.60401522828397e-05];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [1.4878890948170724;1.4392031128470459];
  IW1_1 = [0.072467863579031261 -0.033465039054327989 0.80250838022899018 0.55643451341645822 0.47031727374910498 0.27361609321356045 0.29837707621762399 -0.19837932292577459 -0.47454469505758073 0.3581624476030067 0.21960208645890497 -0.34204236052997761 -0.23175915269427441 -0.5197988140437404 -0.90223202132001945 -1.0530774785860124 -0.93463666154352509 -0.50588624203872168 -0.74043192908467614 -0.11921920556446172 0.12163729923985081 0.42637078218772273 0.71057203593057805 0.51207422044291995 0.46368956587541948 0.28923154552380892 0.29440830211290336 -0.17964016868667024 0.019551030974928406 0.10216715807552042 0.01065938678709619 -0.18030432539240934 0.12715273754188838 -0.15566306158650617 -0.34094531864126038 -0.074217515845973067 -0.02362863875166574 -0.31343188544556927 0.78807952270182013 -0.13576749336664679 0.034230321722001557 0.29176498895099889 -0.0013900099108482293 0.0059065361871274459 0.15421406314503769 -0.29562793334139048 0.25058998571475266 0.13773539550105041 0.20761697057569883 0.05487812777734119 0.2290322183535537 -0.011631473526001322 0.2087864210613973 -0.043942436406814403 0.054212408616148761 0.18644977515537572 -0.024795283747577029 -0.18402275365498219 0.16954481619717049 -0.1383391082739458 0.022955879046180391 -0.069749452562907988 -0.13073793647981732 -0.26453824273302545 -0.013009769935196452 -0.10783741014221257 0.18809553007989424 0.051775088940042117 0.07720041966027319 -0.15080303254780889 0.0072623026707726475 -0.022888177804181627 -0.13509650531064157 0.21118188838040405 -0.12332509704416315 -0.039281218703556949 -0.12547948021624461 -0.11856440084661944 -0.012599612134589608 0.053232782725167418 -0.041891185725510663 -0.070697945076963864 -0.039899756864509854 0.056858966151548784 -0.14248065408144803 0.11527498781643282 -0.023684601154724722 -0.047104960267868839 0.070207965070491651 0.1345449159315559 -0.02075125570403883 0.053328531695836921 0.13186383857609421 -0.0015059484969364649 0.045073380275039467 -0.050945656122786204 -0.027014105116658933 -0.024834084673599002 -0.081410286969283152 -0.013847140471871746 0.14621620469347676 -0.10817541110071141 0.17381131734339253 -0.066184438024457101 0.0090868072983264787 -0.034520579183083294 0.12402660243085242 -0.12315702084305825 -0.046009305952447857 0.011195991110418753 -0.038800272831720876 -0.087993970073904429 0.13779809993623315 0.049152833915679837 0.010670130822078242 0.036182405090420526 -0.096216330917779772 0.07237281204117646 0.0048433418747393488 0.19535625412381111 0.15599824994271955 0.14879711880031327 -0.068036084461149129 0.073583159467551118 0.073820738900818375 0.069019216624192831 0.1264003254865807 0.17611510963165219 -0.011208328240710084 -0.05538318676411829 0.10260030232212781 0.0013827498959071332 0.070695085388703213 0.12939533671042103 0.07673164512952127 -0.036199739318461717 -0.08909455465901743 -0.079840832977124007 0.10180835827282796 -0.0012901394276709829 0.08380959476282189 -0.093618060905002798 -0.083049959478865856 0.15004105151011501 -0.063976596129478494 0.028122261954781444 0.05196382808135807 0.0052219993286825991 0.10244835420496473 -0.051814929297142785 0.0064558520185083543 0.13949532129694164 0.0059757506878363236 0.0096606390948158292 0.023718592780861501 0.0013393560248670916 0.1695877648899923 -0.093756638482385773 -0.023236831910484368 0.0058873592972492105 0.16024955488860865 0.080788506558456458 -0.014524018459135511 0.12656513434890604 0.20627796908162763 0.098276273035541872 0.048544772975178008 0.10186912555161994 -0.083903286182871406 -0.1123192717232358 0.036629128642858381 0.056797375510996791 -0.021437072333076089 0.21325531737032211 0.16210662967874451 -0.028244595224399337 -0.0029732446247071154 -0.13147604768597609 -0.083816821887357285 0.06850600429782383 -0.04742902936519882 -0.02568457167926657 0.21299129925281121 -0.035750014593922533 -0.12295676410060995 0.14672509559573702 -0.088137257308790157 -0.12562666824060423 -0.0008045463599320276 -0.073780844088258662 -0.10303929479262486 0.05928447717922615 -0.022231991445482723 -0.15952906825360777 -0.033436858719116273 0.05188115932441973 0.10669399275707685 0.2149685685226741 -0.076647098632872085 0.039759211390609006 0.017463388331555845 -0.10382237707401869 -0.11010906199975021 0.063239740122320337 0.16709201021089903 0.16265627041771563 0.20416184541489693 -0.097272987793508195 0.049549633368049989 -0.2044736448415215 -0.097474098431709191 0.073294131690631481 0.23030671633379449 -0.22316199305252521 0.012412015486360204 -0.088742415596803362 -0.059471428181312853 0.068124928046498653 0.18338349491299014 -0.23528013535476597 0.15637150658757962 -0.032856039272744948 0.089991754657628514 0.11129320203914764 0.25182660958780773 -0.089368062898280573 0.014831208388001569 -0.017655002927725219 -0.018149361891056745 0.32683651268227409 -0.00064699052066519082 0.048214532844473378 0.86066489723653095 -0.31661195154219174 -0.12318521135243185 0.10974731329262812 -0.40678704145384176 -0.25780393052237544 0.29268557047700094 -0.16101897430991693 -0.077114787349084335 0.039077006025707815 -0.02034951259880452 -0.090208233442323885 0.41416561946417357 0.40815293777179579 0.60394533742719048 0.7708013059324581 0.45232514806919299 0.23408791034825749 0.18011859970743921 0.044079070710913791 -0.60984156966571812 -0.36235101564014599 -0.94181907870158854 -1.2303517433978901 -0.97011275698195454 -0.56865270219854236 -0.17679108279360403 -0.42529548131407285 0.1495371305132332 0.42705886375706009 -0.34206652264544679;-0.38766870977560891 -0.059342643956874469 -0.39597767183535532 -0.124429475663716 -0.32410417860503182 -0.49490626327652248 -0.24054723546992851 -0.28146495908372071 -0.11087898187131999 -0.2446448556272138 -0.024085019948544294 -0.059089104594331671 -0.30857848770544538 -0.065084888516222683 0.015053023003404595 0.10726221811427795 0.16143922003789293 0.52689514159927275 0.40188322000257914 0.29478313561250546 0.30944658546748571 0.52449480388807945 0.40493209951472603 0.31145555028122085 0.47430325774221715 0.71326309267334065 0.81450575102769707 0.39817416295263403 0.53221322107965618 0.43639115504072612 0.18192164462973479 0.47352686066895711 0.40224088933903951 0.21555851019575614 0.045910575782123243 0.10091866943711639 0.17949159426247996 -0.24889530305631788 -0.37620311646879601 -0.41410913662011872 -0.37030421282027409 -0.2424260775768533 -0.35844484169176144 -0.18021722709684521 -0.16381043801057876 -0.26881224508585561 -0.20199123076803269 -0.12627282137329723 -0.33599436648561193 -0.17618954934569131 -0.2928808178647988 0.061028923160451443 -0.19352656207778246 -0.055510452947815793 -0.067833259999704826 -0.055433424218665164 -0.11544132448341247 0.17860159731481495 0.019719887413459952 0.15003335545590862 -0.019851514140179861 0.01713685825691192 -0.0767198820733518 -0.056996075436614821 -0.024203717589807508 0.1136266924245752 0.012476803974117576 -0.058036256102979689 0.01321618713648375 0.080026798246988648 -0.096345775488740737 -0.14100473265898883 0.16676360308444912 0.0085380112890323691 -0.078485193532098332 0.13778539155690003 -0.14419535327947522 0.051145038208731598 0.064751864784805199 0.012517776644641203 0.072845212022058556 0.049212218560225562 -0.042337194224829462 0.097978184983437022 0.1237097282821889 -0.027860287964492272 -0.039384976912472121 -0.016447003740268574 0.086493011871362727 0.031103276313661259 0.091272810475716726 -0.067770459006341738 0.10166633967099131 0.016102523667296514 -0.11822802734725313 -0.031088649363042904 0.080739524890667907 -0.10656791026220563 -0.11019907284961403 -0.076513500991094288 -0.054343148718007758 0.10512995432289667 -0.055797897126143116 -0.098697977571746004 0.0038241714692614105 0.0081249163917004302 -0.11403910738028598 -0.087489066050923422 -0.050860400404933184 -0.16385073856652435 0.073828661239346005 0.10678924619264983 -0.030784848812652744 -0.11007196213547704 0.12453638363447574 -0.028912484949310695 -0.035031135763568338 0.01199333598803479 -0.074887401964026531 -0.18279115610947397 0.10537397668950006 -0.16073682813448784 -0.080316940279400867 -0.15408965161906385 -0.16219215772352691 -0.03833742185363602 0.036510712753188651 0.095548206653013001 0.03366728991001397 -0.18186837611602308 -0.036770021469905406 0.093094405776351127 -0.0040685327042901598 0.067323018377102362 -0.13292751919406326 0.10165049644227346 -0.034068897162860844 0.12309703019772728 -0.041809603179439825 -0.17906029111948921 -0.18413142763947152 -0.023713057554172561 0.072119404313929464 0.041290446668913391 0.018659604602291956 0.11585417369524516 0.061162747318227136 -0.058818701629913983 -0.16143392473726595 -0.13668827279686271 0.078996850984099989 -0.17547050640986617 -0.023542059185850679 0.0048197487604422669 -0.098444654630891734 -0.19175016244117368 -0.13177948153311139 -0.066556316027508708 -0.030437478843437185 0.11412406951107855 -0.16716386895763752 -0.10472354557639747 -0.011023068889446049 -0.048500595659109282 -0.088370099578483741 0.028572911891204752 0.037061449322542635 -0.12443913033543116 0.12461068106949799 -0.17577667511208495 0.091107872110917379 -0.11537341040045802 -0.035144562293958653 -0.12895491298539236 0.016478857665900612 -0.12125894501069967 -0.02989974726602591 0.0054857306758247983 -0.056558149564796351 -0.099626691710296644 0.028639938711540108 0.089088688835664009 0.058353753119832721 0.0078863430531807936 0.10859896031650429 0.078724008284130467 0.015817705643679618 -0.01633474181589608 -0.10899884059198869 -0.035770087595974845 0.025773316531517446 0.06048212816953738 0.061784226299024241 -0.13573615029719988 0.012907413485054496 0.0045989199207413199 0.077245228195555057 -0.0055912691878717558 -0.10259517257197479 0.078736594247361832 -0.14493087913349859 -0.12176019457473405 0.036028771361197046 -0.019486811007380898 0.021232355641477497 0.028284584580643854 -0.081062610960143375 -0.023131682617910102 0.069879277545683308 -0.030208362686036397 0.08219573541836836 -0.013744670972737565 0.060126491249221023 0.092578474824124915 0.13267294133910854 0.034930891628385807 -0.058332346003472924 -0.035350473923570054 -0.068923519381526424 -0.14553499327038538 -0.14189312157658279 -0.10068090280374584 -0.29942069321876824 -0.22143249802492376 -0.16441478495837661 -0.2817951709514378 -0.27955965698508334 -0.18210155025219163 -0.30603085870483526 -0.087214177285703368 -0.26578007158740363 -0.39905305709192174 -0.30396099900274171 -0.23650322268507226 0.026432097047581097 0.17804269464934627 0.23467437198757179 0.36417328748286187 0.24097780237176697 0.36566228236060733 0.27359820083876474 0.42557539773635861 0.63648395954998094 0.48926816465087669 0.6906671587265456 0.62618459898612733 0.37182396180225652 0.33379933500277115 0.26403502035892462 0.36788063911642355 0.26431601488504797 0.25505319119411757 0.51381877321428804 0.27966218029531814 0.34734118358242699 0.06549605830340989 0.20348171557457878 -0.16681915772502784 -0.21336119673602241 -0.007292624854246978 -0.081866449111171552 -0.37162097707535519 -0.33869989105278636];
  
  % Layer 2
  b2 = [-0.96044804927332827;-1.3417410336084954;1.383126401207742;-0.99090647026071355];
  LW2_1 = [3.1398977993103712 3.8648305418292974;-3.8393889798878886 -2.5264070197264368;-1.5341046812406756 1.5185944740089021;2.685129227468412 -2.4080656241830032];
  
  % ===== SIMULATION ========
  
  % Format Input Arguments
  isCellX = iscell(X);
  if ~isCellX, X = {X}; end;
  
  % Dimensions
  TS = size(X,2); % timesteps
  if ~isempty(X)
    Q = size(X{1},2); % samples/series
  else
    Q = 0;
  end
  
  % Allocate Outputs
  Y = cell(1,TS);
  
  % Time loop
  for ts=1:TS
  
    % Input 1
    Xp1 = mapminmax_apply(X{1,ts},x1_step1_gain,x1_step1_xoffset,x1_step1_ymin);
    
    % Layer 1
    a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*Xp1);
    
    % Layer 2
    a2 = softmax_apply(repmat(b2,1,Q) + LW2_1*a1);
    
    % Output 1
    Y{1,ts} = a2;
  end
  
  % Final Delay States
  Xf = cell(1,0);
  Af = cell(2,0);
  
  % Format Output Arguments
  if ~isCellX, Y = cell2mat(Y); end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings_gain,settings_xoffset,settings_ymin)
  y = bsxfun(@minus,x,settings_xoffset);
  y = bsxfun(@times,y,settings_gain);
  y = bsxfun(@plus,y,settings_ymin);

% Competitive Soft Transfer Function
function a = softmax_apply(n)
  nmax = max(n,[],1);
  n = bsxfun(@minus,n,nmax);
  numer = exp(n);
  denom = sum(numer,1); 
  denom(denom == 0) = 1;
  a = bsxfun(@rdivide,numer,denom);

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n)
  a = 2 ./ (1 + exp(-2*n)) - 1;
