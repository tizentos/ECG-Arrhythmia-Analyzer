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
end
