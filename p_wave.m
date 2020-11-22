function [pr_int,pr_interval,p_peak,p_sample,p_int,locs_Pwave]=p_wave(ECG_data,fs,locs_Qwave)
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
        locs_Pwave(i)=ind1(i)+origin_zero1(i);
    elseif ind==2
        pr_int=P2_int;
        pr_interval=(length(pr_int(i,:))/fs)*1000;
        p_sample(i,:)=ECG_data(origin_zero2(i):round(length(pr_int(i,:))/6):pr_int(i,end));
        p_int(i,:)=origin_zero2(i):round(length(pr_int(i,:))/6):pr_int(i,end);
        locs_Pwave(i)=ind2(i)+origin_zero2(i);
    end
end
% p_peak=sum(p_peak)/length(p_peak);