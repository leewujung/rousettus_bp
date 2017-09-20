clear;

date_data='20150824';
trl_num='01';

fn=['..\Test\' date_data '\Session 1\Trial' trl_num '.c3d'];
[point_array, frame_rate, trig_rt, trig_sig, start_f, end_f] = lc3d( fn );
point_names=cellfun(@(c) c.name,point_array,'uniformoutput',0);

% if not using calib mics...
if datenum(date_data,'yyyymmdd') < 736200 %20150824
  MN=load('mic_names.mat');
  no_calib=1;
else % using calib mics:
  MN=load('mic_names_with_calib.mat');
  no_calib=0;
end
mic_names=MN.names_mics;

mic_pos=nan(length(mic_names),3);
for mm=1:length(mic_names)
  imic=find(strcmp(point_names,mic_names{mm}));
  if ~isempty(imic)
    pmic=point_array{imic}.traj./1e3;
    indx=find(pmic(:,1),1);
    mic_pos(mm,:)=pmic(indx,:);
  end
end

figure(2); clf; set(gcf,'pos',[10 300 912 700]), hold on;
cols='rgb';
if no_calib
  for mm=1:3
    scatter3(mic_pos(mm:3:end,3),mic_pos(mm:3:end,1),mic_pos(mm:3:end,2),...
      cols(mm))
  end
else
  for mm=1:3
    scatter3(mic_pos(mm:3:end-4,3),mic_pos(mm:3:end-4,1),mic_pos(mm:3:end-4,2),...
      cols(mm))
  end
  scatter3(mic_pos(end-3:end,3),mic_pos(end-3:end,1),mic_pos(end-3:end,2),'k')
  tip_indx=[(1:3:32*3)+2 32*3+2 32*3+4];
  tips=mic_pos(tip_indx,:);
  text(tips(:,3),tips(:,1),tips(:,2),num2str((1:length(tips(:,1)))'))
end
view(3)
axis equal, grid on  
view(-57,22)
% scatter3(mic_pos(:,1),mic_pos(:,2),mic_pos(:,3));

save(['..\mic_pos\mic_pos_' date_data '_S1_' trl_num '.mat'],'mic_pos')

% names_mics=cellfun(@(c) c.name,mics,'uniformoutput',0);
% save('mic_names_with_calib.mat','names_mics');