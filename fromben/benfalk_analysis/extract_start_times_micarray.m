clear
dirname='..\mic_recordings\20151218_rousettus_matfile\';
files=dir([dirname '*.mat']);

ST=cell(length(files),2);
for ff=1:length(files)
  m=matfile([dirname files(ff).name],'writable',false);
  ST{ff,1}=files(ff).name;
  ST{ff,2}=m.start_time;
end

xlswrite([dirname 'file_start_times.xlsx'],ST)


