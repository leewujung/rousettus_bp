%copy files into new file name format and into single directory
clear;
new_top_dir='..\mic_data\';
top_dir='..\mic_recordings\';

type_file='mic_data';

if strcmp(type_file,'mic_data_detect')
  ideal_split_len=6;
else
  ideal_split_len=5;
end

DD=dir(top_dir);
DD=DD(3:end);
DD=DD([DD.isdir]);
for dd=1:length(DD)
  BD=dir([top_dir DD(dd).name]);
  BD=BD(3:end);
  BD=BD([BD.isdir]);
  for bd=1:length(BD)
    FD=dir([top_dir DD(dd).name '\' BD(bd).name '\*' type_file '.mat']);
    for fd=1:length(FD)
      orig_fn=[top_dir DD(dd).name '\' BD(bd).name '\' FD(fd).name];
      [path,fn]=fileparts(orig_fn);
      D=strsplit(path,'\');
      datename=D{end-1};
      E=strsplit(D{end},'_');
      species=E{1};
      if strcmp(species,'rouesttus')
        species='rousettus';
      elseif strcmp(species,'epteiscus')
        species='eptesicus';
      end
      band=E{2};
      
      C=strsplit(fn,'_');
      if length(C)==ideal_split_len
        trialnum=C{3};
      else
        bandtrial=C{2};
        if isempty(find(isstrprop(C{2}, 'digit')==0,1)) && isempty(find(strfind(C{2},band),1))
          trialnum=C{2};
        else
          F=strsplit(C{2},band);
          trialnum=F{2};
        end
      end
      trialnum=num2str(str2double(trialnum),'%2.2d');
      exportfn=[new_top_dir ...
        strjoin({species,datename,band,trialnum,type_file},'_') '.mat'];
      if ~exist(exportfn,'file')
        copyfile(orig_fn,exportfn)
      end
    end
  end
end
