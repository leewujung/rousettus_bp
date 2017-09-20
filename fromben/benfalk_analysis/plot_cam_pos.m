figure(1); 
plot_cam_nums=1;
plot_cam_dirs=1;

pname='F:\small_space_beampattern\Test\20150825\Session 1\';
fname='Trial02.xcp';

[cam_p,cam_o,cam_n]=extract_cam_positions([pname fname]);

cam_p=cam_p./1e3;

hold on;
plot3(cam_p(:,3),cam_p(:,1),cam_p(:,2),...
  'sk','markerfacecolor','r')
grid on;
axis equal

if plot_cam_dirs
  %plot direction vector from xyz position of camera
  A=[0 0 1];
  B=[0 1 0];
  C=[1 0 0];
  hold on;
  plot3([0 A(3)],[0 A(1)],[0 A(2)],'g')
  plot3([0 B(3)],[0 B(1)],[0 B(2)],'r')
  plot3([0 C(3)],[0 C(1)],[0 C(2)],'b')
  for k=1:size(cam_p,1)
    c = cam_p(k,:);
    q=cam_o(k,:);
    
%     R = quat2rotm(q([4 1:3]));
    R = quat2rot(q([4 1:3]));
%     eul = rotm2eul(R);
    
%     ax_of_R=[R(3,2)-R(2,3) R(1,3)-R(3,1) R(2,1)-R(1,2)];
    
%     E=(ax_of_R-c).*.2;
%     E=A*R.*.2;
    E=R(3,:)*.4;
    quiver3(c(3),c(1),c(2),E(3),E(1),E(2),0,'MaxHeadSize',.5,'linewidth',2,'color','r')
    
%     plot3([c(3) E(3)+c(3)],...
%       [c(1) E(1)+c(1)],...
%       [c(2) E(2)+c(2)],'k')
  end
end

if plot_cam_nums
  text(cam_p(:,3),cam_p(:,1),cam_p(:,2),...
    num2str(cam_n))
end


view(3); axis equal; grid on; 