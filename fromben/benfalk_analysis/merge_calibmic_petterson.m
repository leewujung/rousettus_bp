%merge_calibmic_petterson
clear
bat='36134'; %(session 2)

dates={'20150825','20150910'};
trials={'2','3'};

for dd=1:length(dates)
  for tt=1:length(trials)
    pet_apath=['..\mic_recordings\' dates{dd} '\rousettus_' bat '_matfile\'];
    
    cal_apath=['..\mic_recordings\' dates{dd} '_calib_mic\rousettus_' bat '_matfile\'];
    
    if isequal(dates{dd},'20150825')
      insert='_';
    else
      insert='';
    end
    
    afname=['rousettus_' bat '_' trials{tt} '_mic_data_detect.mat'];
    pet_data=load([pet_apath afname]);
    
    if pet_data.num_ch_in_file < 34
      calafname=['rousettus_' bat insert trials{tt} '_mic_data.mat'];
      cal_data=load([cal_apath calafname]);
      
      y=resample(cal_data.sig,cal_data.fs/4,cal_data.fs);
      y(end+1,:)=nan;
      pet_data.sig(:,end+1:end+2)=y;

      sig=pet_data.sig;
      num_ch_in_file=34;

%       save([pet_apath afname],'sig','num_ch_in_file','-append')
    end
    
%     yy=[];
%     for ch=1:size(cal_data.sig,2)
%       yy(:,ch)=cal_data.sig(3:4:end,ch);
%     end
%     yyy=downsample(cal_data.sig,4);
    
    
%     figure(1), clf;
%     plot((1:length(cal_data.sig))/cal_data.fs,cal_data.sig)
%     hold on;
%     
%     plot((1:length(y))/(cal_data.fs/4),y)
%     plot((1:length(yy))/(cal_data.fs/4),yy)
%     plot((1:length(yyy))/(cal_data.fs/4),yyy)
%     legend(num2str((1:6)'))
  end
  
end