function [qt_int,t_sample,t_peak,qt_interval,rt_interval,sample_int,locs_Twave]=t_wave(ECG_data,fs,locs_Qwave,locs_Rwave)
%  T-onset and off-set must be considered
%  the peak between duration=T_onset-T_offset is the T-wave
Q_T1=ceil(0.44*fs);
Q_T2=ceil(0.4*fs);
Q_Tm=ceil(0.36*fs);
% intervals
if (locs_Qwave(end)+Q_T1>ECG_data(end))
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
        locs_Twave(i)=(int1(i)-T_int)+ind1(i);
    elseif ind==3
        qt_int(i,:)=(int2(i)-T_int):int2(i);
        qt_interval=(Q_T2/fs)*1000;
        RT_int(i)=int2(i)-locs_Rwave(i);
        rt_interval(i)=(RT_int(i)/fs)*1000;
        t_sample(i,:)=ECG_data(int2(i)-T_int:round(length(Q_section2)/6):int2(i));
        sample_int(i,:)=int2(i)-T_int:round(length(Q_section2)/6):int2(i);
        locs_Twave(i)=(int2(i)-T_int)+ind2(i);
    elseif ind==1
        qt_int(i,:)=(intm(i)-T_int):intm(i);
        qt_interval=(Q_Tm/fs)*1000;
        RT_int(i)=intm(i)-locs_Rwave(i);
        rt_interval(i)=(RT_int(i)/fs)*1000;
        t_sample(i,:)=ECG_data(intm(i)-T_int:round(length(Q_sectionm)/6):intm(i));
        sample_int(i,:)=intm(i)-T_int:round(length(Q_sectionm)/6):intm(i);
        locs_Twave(i)=(intm(i)-T_int)+indm(i);
    end    
end
% t_peak=sum(t_peaks)/length(t_peaks);
rt_interval=sum(rt_interval)/length(rt_interval);
%what about average of qt_interval
