% assign_mics
%takes the mic order as marked in vicon and matches to channel layout
%this is just a plotting and helper function, you have to do the assignment
%manually

clear;

date_data='20150824';
trl_num='01';

% if not using calib mics...
if datenum(date_data,'yyyymmdd') < 736200 %20150824
  no_calib=1;
  tip_indx=[(1:3:32*3)+2];
else % using calib mics:
  no_calib=0;
  tip_indx=[(1:3:32*3)+2 32*3+2 32*3+4];
end

load(['..\mic_pos\mic_pos_' date_data '_S1_' trl_num '.mat'])
tips=mic_pos(tip_indx,:);

figure(1); set(gcf,'pos',[10 300 912 700]), clf; hold on;
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
end
scatter3(tips(:,3),tips(:,1),tips(:,2),'.k')
text(tips(:,3),tips(:,1),tips(:,2),num2str((1:length(tips(:,1)))'))
axis equal, grid on  
view(-69.2,22)

if exist(['..\mic_pos\mic_order_' date_data '_S1_' trl_num '.mat'],'file')
  load(['..\mic_pos\mic_order_' date_data '_S1_' trl_num '.mat'])
else
  if no_calib
    mic_order=nan(32,1);
  else
    mic_order=nan(32+2,1);
  end
end

disp('enter your mic order')
pause;

save(['..\mic_pos\mic_order_' date_data '_S1_' trl_num '.mat'],'mic_order')
