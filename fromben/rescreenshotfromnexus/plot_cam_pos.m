plot_cam_nums=1;
plot_cam_dirs=1;

pname='D:\small_space_beampattern\Test\20150825\Session 1\';
fname='Trial02.xcp';

[cam_p,cam_o,cam_n]=extract_cam_positions([pname fname]);

cam_p=cam_p./1e3;

hold on;
plot3(cam_p(:,3),cam_p(:,1),cam_p(:,2),...
  'sk','markerfacecolor','r')
grid on;
axis equal

if plot_cam_dirs
  %convert quaternion to direction vector
  %plot direction vector from xyz position of camera
  
  %cam_o
end

if plot_cam_nums
  text(cam_p(:,3),cam_p(:,1),cam_p(:,2),...
    num2str(cam_n))
end

