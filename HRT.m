function heart_rate=HRT(locs_Rwave,f_sampling)
avg=[];
j=1;
for i=2:length(locs_Rwave)
    avg(j)=locs_Rwave(i)-locs_Rwave(i-1);
    j=j+1;
end
heart_rate=(sum(avg)/length(avg))*(1/f_sampling);
end

