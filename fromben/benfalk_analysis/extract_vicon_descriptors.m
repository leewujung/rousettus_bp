%extract the description text from each of the calibration files and store
%them in a central location

descfileID = fopen('calib_descriptors.txt','w');

top_dir='..\Test\';
DD=dir(top_dir);
DD=DD(3:end);
DD=DD([DD.isdir]);
for dd=1:length(DD)
  BD=dir([top_dir DD(dd).name]);
  BD=BD(3:end);
  BD=BD([BD.isdir]);
  if ~isempty(BD)
    for bd=1%:length(BD)
      FD=dir([top_dir DD(dd).name '\' BD(bd).name '\Trial*.enf']);
      for fd=1:length(FD)
        enf_file=[top_dir DD(dd).name '\' BD(bd).name '\' FD(fd).name];
        
        fileID = fopen(enf_file,'r');
        E=textscan(fileID,'%s','Delimiter','\n');
        fclose(fileID);
        
        rownum= ~cellfun('isempty',strfind(E{1},'DESCRIPTION'));
        F=strsplit(E{1}{rownum},...
          'DESCRIPTION=');
        if length(F)>1
          desc=F{2};
        else
          desc='';
        end
        
        fprintf(descfileID,'%s: %s\n',[DD(dd).name '\' FD(fd).name],desc);
      end
    end
  end
end

fclose(descfileID);
