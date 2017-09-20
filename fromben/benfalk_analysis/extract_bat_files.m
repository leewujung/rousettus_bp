%extract bat data across files
clear;
top_dir='..\Test\';

DD=dir(top_dir);
DD=DD(3:end);
DD=DD([DD.isdir]);
for dd=1:length(DD)
  BD=dir([top_dir DD(dd).name]);
  BD=BD(3:end);
  BD=BD([BD.isdir]);
  for bd=1:length(BD)
    FD=dir([top_dir DD(dd).name '\' BD(bd).name '\*.c3d']);
    for fd=1:length(FD)
      orig_fn=[top_dir DD(dd).name '\' BD(bd).name '\' FD(fd).name];
      
      extract_bat(orig_fn);
    end
  end
end
