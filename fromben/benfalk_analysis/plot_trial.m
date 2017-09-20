
clear;

date_data='20150825';
trl_num='02&04'; %for mics
% date_data='20150910';
% trl_num='04'; %for mics


%load mics
load(['..\mic_pos\mics_ch_order_' date_data '_S1_' trl_num '.mat']);
load(['..\mic_pos\mic_pos_' date_data '_S1_' trl_num '.mat']);

%load bat
load('..\\bat_pos\20150825_S2_02.mat');
% load('..\\bat_pos\20150910_S2_02.mat');
bat=bat_pos{1};
frames=find(isfinite(bat_pos{1}(:,1)));





v=VideoWriter('bat2mics_trial.mp4','mpeg-4');
v.FrameRate=15;
v.Quality=90;
open(v);
figure(1), clf; set(gcf,'pos',[10 40 800 600])
for ff=frames(180:end)'
  cla;
  plot3(bat(:,3),bat(:,1),bat(:,2),'-','color',[.5 .5 .5])
  plot3(bat(ff,3),bat(ff,1),bat(ff,2),'ok');
  hold on;
  plot3(tips(:,3),tips(:,1),tips(:,2),'ok')
  plot3([zeros(length(mic_vec),1)'+bat(ff,3); tips(:,3)'],...
    [zeros(length(mic_vec),1)'+bat(ff,1); tips(:,1)'],...
    [zeros(length(mic_vec),1)'+bat(ff,2); tips(:,2)'],'k')
  grid on; axis equal
  view(-57,21)
  
  writeVideo(v,getframe(gcf));
end
close(v)



cols='rgb';
figure(1), clf; set(gcf,'pos',[10 40 800 600]); hold on;
for mm=1:length(markers)
  plot3(bat_pos{mm}(:,3),bat_pos{mm}(:,1),bat_pos{mm}(:,2),[cols(mm) '-'],'linewidth',2)
end
plot3(tips(:,3),tips(:,1),tips(:,2),'ok')
for mm=1:3
  scatter3(mic_pos(mm:3:end-4,3),mic_pos(mm:3:end-4,1),mic_pos(mm:3:end-4,2),...
    cols(mm))
end
scatter3(mic_pos(end-3:end,3),mic_pos(end-3:end,1),mic_pos(end-3:end,2),'k')

axis equal, grid on, view(3)



v=VideoWriter(['bat2mics_angles_' date_data '.mp4'],'mpeg-4');
v.FrameRate=15;
v.Quality=90;
open(v);
set(gcf,'pos',[10 40 800 600])
for ff=frames(100:end-50)'
  cla;
  bat_mic_vec=tips(:,[3 1 2])-repmat(bat(ff,[3 1 2]),34,1);
  [az,el,R] = cart2sph(bat_mic_vec(:,1),bat_mic_vec(:,2),bat_mic_vec(:,3));

  scatter(az*180/pi,el*180/pi,'+k')
  axis([-90 90 -90 90])
  grid on;
  set(gca,'fontsize',16)
  
  writeVideo(v,getframe(gcf));
end
close(v)