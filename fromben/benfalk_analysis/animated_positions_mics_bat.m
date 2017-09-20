clear

load('D:\small_space_beampattern\mic_pos\mic_pos_20150824_S1_01.mat')

colordef black
figure(2); clf; set(gcf,'pos',[10 300 912 700],'color','black'), hold on;
% cols='rgb';
cc=colormap;

for mm=1:3
  scatter3(mic_pos(mm:3:end-4,3),mic_pos(mm:3:end-4,1),mic_pos(mm:3:end-4,2),...
    [],'filled')
end
scatter3(mic_pos(end-3:end,3),mic_pos(end-3:end,1),mic_pos(end-3:end,2),'k','filled')
tip_indx=[(1:3:32*3)+2 32*3+2 32*3+4];
tips=mic_pos(tip_indx,:);
% text(tips(:,3),tips(:,1),tips(:,2),num2str((1:length(tips(:,1)))'))

view(3)
axis equal, grid on
view(-60,12)


files=dir('..\bat_pos\eptesicus_20150824*.mat');
for ff=1:length(files)
  load(['..\bat_pos\' files(ff).name]);
  bat=bat_pos{1};
  plot3(bat(:,3),bat(:,1),bat(:,2),'color',cc(end-10,:))
end

v = VideoWriter('F:\newfile.avi','uncompressed avi');
v.FrameRate=20;
open(v)
for k=[-60:60 60:-1:-60]
  view(k,12)
%   drawnow
  writeVideo(v,getframe(gcf))
end
close(v)
compress_video([v.Path v.Filename],1,0)