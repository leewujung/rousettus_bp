clear;
useline=1;
DIAG=0;
bat_pos_dir='..\bat_pos\';
files=dir([bat_pos_dir '*.mat']);
sideLeft=nan(length(files),820);
for k=1:length(files)
  load([bat_pos_dir files(k).name])
  
  frames_w_alldata=find(isfinite(bat_pos{1}(:,1)) & isfinite(bat_pos{2}(:,1)));
  
  sideLeft(k,:)=determine_side(bat_pos{1}(:,[3,1,2]),bat_pos{2}(:,[3,1,2]),...
    frames_w_alldata,useline,DIAG);
  if DIAG
    title(files(k).name,'Interpreter','none')
  end
end

wrong_markings = find(mode(sideLeft,2)==-1);
% wrong_markings2 = find(nanmean(sideLeft,2)<0);

wrong_marked_fnames={files(wrong_markings).name}'