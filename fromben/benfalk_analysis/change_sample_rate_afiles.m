clear;
mat_dir='G:\small_space_beampattern\mic_recordings\20150910\rousettus_39184_matfile\';
% mat_dir='G:\small_space_beampattern\mic_recordings\20150910\rousettus_36134_matfile\';

files=dir([mat_dir '*.mat']);
for f=1:length(files)
  m=matfile([mat_dir files(f).name],'Writable',true);
  m.fs=250e3;
end

clear;