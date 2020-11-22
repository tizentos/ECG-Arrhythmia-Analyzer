function [class_vect,beat]=analyzeECG(ecg_signal)
sig=ecg_signal;
n=ceil(log2(length(sig)));
t=10;    %duration of signal
%%f=length(sig)/t;   %%sampling frequency that samples per second
f=360;
freq_scaling=f*(0:(2^n)-1)/2^n;
%%detrending
[p,s,mu] = polyfit((1:numel(sig)),sig,6);
f_y = polyval(p,(1:numel(sig)),[],mu);
% h=fir1(1000,4/(f*2),'high');
% l=fir1(1000,100/(f*2),'low');
% filt_ECG=filter(h,1,sig);
% filt_ECG=filter(l,1,filt_ECG);
ECG_data = sig-f_y;
% end detrending
% plot(ECG_data)
%%filtering
%  ECG_data=medfilt1(ECG_data,10);
% ECG_bpf_struct=load('ECG_bpf.mat');
% ECG_bpf=ECG_bpf_struct.ECG_bpf;
% ECG_data=filter(ECG_bpf,ECG_data);
for i=1:10
    ECG_data=medfilt1(ECG_data,4);
end
% filtering end
%%analysis proper
f_comp=fft(ECG_data,2^n);
sig_energy=abs(f_comp.*f_comp);
figure('Name','ECG ANALYSIS')
subplot(2,2,3)
plot(freq_scaling,abs(f_comp))
title('signal fft')
grid on
xlabel('frequency(Hz)'); ylabel('voltage level')
subplot(2,2,4)
plot(freq_scaling,sig_energy)
title('signal energy')
grid on
xlabel('frequency(Hz)'); ylabel('voltage level')
%% my algorithm on peak finding
[locs_Qwave,locs_Rwave,locs_Swave,qrs_area,qrs_duration,q_samples]...
    =qrs_comp(ECG_data,f);
%% heart rate calculation and QRS complex extraction
 heart_rate=HRT(locs_Rwave,f);
 R_R=heart_rate*1000;
 bpm=(1/heart_rate)*60;
fprintf('\tR-R : %.0fms \n\tHeart Rate :%.0fbpm\n\tQRS : %.0fms\n',R_R,bpm,round(qrs_duration))  
%END of QRS evaluation
%% T-wave detection
[~,~,t_peak,qt_interval,rt_interval,sample_int,locs_Twave]=t_wave(ECG_data,f,locs_Qwave,locs_Rwave);
%[~,t_peak,qt_interval,rt_interval,sample_int]=t_wave2(ECG_data,f,locs_Qwave,locs_Rwave);
%% P-WAVE DETECTION
 [~,pr_interval,p_peak,~,p_int,locs_Pwave]=p_wave(ECG_data,f,locs_Qwave);
%[~,pr_interval,p_peak,~,p_int,locs_Pwave]=p_wave2(ECG_data,f,locs_Qwave);
%%plotting
subplot(2,2,[1 2])
 hold on
 plot(ECG_data)
 plot(locs_Swave,ECG_data(locs_Swave),'rs','MarkerFaceColor','b');
 plot(locs_Rwave,ECG_data(locs_Rwave),'rs','MarkerFaceColor','g');
 plot(locs_Qwave,ECG_data(locs_Qwave),'rs','MarkerFaceColor','r'); 
grid on
xlabel('sample number'); ylabel('voltage level')
 figure
hold on
plot(ECG_data)
plot(locs_Pwave,ECG_data(locs_Pwave),'rs','MarkerFaceColor','g');
plot(locs_Twave,ECG_data(locs_Twave),'rs','MarkerFaceColor','b');
plot(locs_Twave,t_peak,'rs','MarkerFaceColor','b');
plot(locs_Pwave,p_peak,'rs','MarkerFaceColor','r');  % CHANGE
%plot(q_sample_int,q_samples,'rs','MarkerFaceColor','g');
%RT_INT TO INTERVAL
grid on
xlabel('sample number'); ylabel('voltage level')
%%beat Extraction
%%beat=beatExtractor(locs_Rwave,p_int,sample_int,ECG_data);
%%classification vector
class_vect=[qrs_duration,qrs_area,pr_interval,rt_interval,qt_interval,t_peak,p_peak,sig_energy];
beat=beatExtractor(locs_Rwave,locs_Pwave,locs_Twave,ECG_data);