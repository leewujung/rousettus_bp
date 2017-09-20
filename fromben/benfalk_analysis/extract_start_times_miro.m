clear
dirname='F:\Miro_20151218\';
files=dir([dirname '*.xml']);

ST=cell(length(files),2);
for ff=1:length(files)
  s = xml2struct( [dirname files(ff).name] );
  start_time = [s.chd.CineFileHeader.TriggerTime.Date.Text ' ' s.chd.CineFileHeader.TriggerTime.Time.Text];
  
  ST{ff,1}=files(ff).name;
  ST{ff,2}=start_time;
end

xlswrite([dirname 'file_start_times.xlsx'],ST)