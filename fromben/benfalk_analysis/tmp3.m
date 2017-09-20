clear
D = '..\bat_pos\';
ff = dir([D,'*_bat_pos.mat']);
for iF=1:length(ff)
  fname = ff(iF).name;
  fname_post = strsplit(fname,'_');
  nname = strjoin(fname_post([1 3 2 4:end]),'_');
  nname = regexprep(nname,'__','_');
  movefile([D fname],[D nname])
end