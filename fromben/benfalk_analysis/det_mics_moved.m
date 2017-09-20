%how much do microphones move?
clear
mic_pos_dir = '..\mic_pos\';

date='20150819';

basefn=['mic_pos_' date '_S1_02.mat'];
basemics = load([mic_pos_dir basefn]);
BMP=basemics.mic_pos;

files=dir([mic_pos_dir 'mic_pos_' date '_S1_*.mat']);
findx=~ismember({files.name},basefn);

MD=[];

for fi=find(findx)
  load([mic_pos_dir files(fi).name])
  
  for mp=find(isfinite(mic_pos(:,1)'))
    MP=mic_pos(mp,:);
    
    DD=distance(BMP,MP);
    [MD(fi,mp),id]=min(DD);
    
    
  end
end