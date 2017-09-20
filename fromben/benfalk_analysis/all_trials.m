clear; 

vicon_dir='..\Test\';
dirs=dir(vicon_dir);
dirs=dirs(3:end);
dirs=dirs([dirs.isdir]);

for dd=1:length(dirs)
  subdir=dirs(dd).name;
  dirs=dir([vicon_dir subdir]);
  dirs=dirs(3:end);
  dirs=dirs([dirs.isdir]);
  