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
qrs_duration=((locs_Swave-locs_Qwave)/fs)*1000;  % convert to seconds using 1000
qrs_duration=sum(qrs_duration)/length(qrs_duration);
qrs_area=0.5.*R_amp.*(qrs_duration/1000);
qrs_area=sum(qrs_area)/length(qrs_area);
%% for q_samples
sample_space=locs_Swave-locs_Qwave;
q_samples=ECG_data(locs_Swave:round(sample_space/6):locs_Qwave);
end