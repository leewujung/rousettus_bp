clear
top_dir='E:\0_WJLEE_DATA\small_space_beampattern\mic_data\';
fd=dir([top_dir 'eptesicus_20150824_calib_mic_LB53_*_mic_data.mat']);

for ff=1:length(fd)
  disp(['on file: ' regexprep(fd(ff).name,'calib_mic_','')])
  pet_fn=[top_dir regexprep(fd(ff).name,'calib_mic_','')];
  detect_fn=[pet_fn(1:end-4) '_detect.mat'];
  if exist(pet_fn,'file')
    try
      pet_data=load([top_dir pet_fn]);
    catch
      disp(['File cannot be read: ' pet_fn])
      n=0;
      continue
    end
    
    if ~isempty(pet_data.sig)
      cal_data=load([top_dir fd(ff).name]);
      y=cal_data.sig;
      if length(y) == length(pet_data.sig)
        pet_data.sig(:,end+1:end+size(y,2))=y;

        sig=pet_data.sig;
        save(pet_fn,'sig','-append')
        delete([top_dir fd(ff).name]);

        if exist(detect_fn,'file')
          num_ch_in_file=size(pet_data.sig,2);
          save(detect_fn,'num_ch_in_file','-append')
        end
      end
    end
  end
end