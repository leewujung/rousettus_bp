%remove sig from data_detect files
clear;
top_dir='..\mic_recordings\';

DD=dir(top_dir);
DD=DD(3:end);
DD=DD([DD.isdir]);
for dd=1:length(DD)
  BD=dir([top_dir DD(dd).name]);
  BD=BD(3:end);
  BD=BD([BD.isdir]);
  for bd=1:length(BD)
    FD=dir([top_dir DD(dd).name '\' BD(bd).name '\*data_detect.mat']);
    for fd=1:length(FD)
      sig_fname=[top_dir DD(dd).name '\' BD(bd).name '\' FD(fd).name(1:end-11) '.mat'];
      if exist(sig_fname,'file')
        msig=matfile([top_dir DD(dd).name '\' BD(bd).name '\' FD(fd).name]);
        Ssig=whos(msig);
        if ~isempty(find(strcmp({Ssig.name},'sig'), 1))
          m=matfile([top_dir DD(dd).name '\' BD(bd).name '\' FD(fd).name],...
            'Writable',true);
          S=whos(m);
          if ~isequal(S(strcmp({S.name},'sig')).size,[0 0]) || ...
              ~isequal(S(strcmp({S.name},'sig_t')).size,[0 0])
            m.sig=[];
            m.sig_t=[];
          end
        end
      end
    end
  end
end
