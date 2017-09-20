%merge sig from data_detect files
clear;
top_dir='..\mic_data\';

type_file='calib_mic';

FD=dir([top_dir '*' type_file '*.mat']);

tic; to=[]; n=0;
for fd=1:length(FD)
  %displaying progress
  if mod(fd,round(length(FD)/100))==0
    to(end+1)=toc;
    fprintf(repmat('\b',1,n));
    sec_left=mean(to)*( (length(FD)-fd)/length(FD)*100);
    min_left=floor(sec_left/60);
    sec_left=round(rem(sec_left,60));
    msg=[num2str(fd/length(FD)*100,'%2.0f') '% done, time left: ' ...
      num2str( min_left ,'%2.0f' ) ' m ' num2str( sec_left ,'%2.0f' ) ' s' ];
    fprintf('%s',msg);
    n=numel(msg);
    tic
  end
  
  cal_fn=[top_dir FD(fd).name];
  pet_fn=[top_dir regexprep(FD(fd).name,'calib_mic_','')];
  detect_fn=[pet_fn(1:end-4) '_detect.mat'];
  if exist(pet_fn,'file')
    try
      pet_data=load(pet_fn);
    catch
      disp(['File cannot be read: ' pet_fn])
      n=0;
      continue
    end
  else
    disp(['File does not exist: ' pet_fn])
    continue;
  end

  if ~isempty(pet_data.sig)
    cal_data=load(cal_fn);
    y=resample(cal_data.sig,cal_data.fs/4,cal_data.fs);
    y=[zeros(1,size(y,2)); y];

    if length(y) == length(pet_data.sig)
      pet_data.sig(:,end+1:end+size(y,2))=y;

      sig=pet_data.sig;
      save(pet_fn,'sig','-append')
      delete(cal_fn);

      if exist(detect_fn,'file')
        num_ch_in_file=size(pet_data.sig,2);
        save(detect_fn,'num_ch_in_file','-append')
      end
    end
  end
end
