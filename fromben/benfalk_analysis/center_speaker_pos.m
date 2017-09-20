%calc_circle_centers
clear;
close all;
DIAG=1;

data_date='20150826';

load(['speaker_pos_' data_date '.mat']);
D=[]; XYZc=[];
for sp=1:length(speaker_pos)
  XYZ=speaker_pos{sp};
  scentroid=mean(XYZ,1);
  
  XYZ0 = bsxfun(@minus,XYZ,scentroid);
  [U,S,V] = svd(XYZ0,0);
  P = V(:,3);
  
  v=[0 0 1];
  %rotating along one axis
  th=atan2(norm(cross(P,v)),dot(P,v));
  rot_mat=[cos(th) 0 sin(th); 0 1 0;-sin(th) 0 cos(th)]; % rot_mat_y
  for k=1:length(XYZ0)
    XYZr(k,:)=rot_mat*(XYZ0(k,:))';
  end
  [xc,yc,R,a]=circfit(XYZr(:,1),XYZr(:,2));
%   par=CircleFitByTaubin([XYZr(:,1),XYZr(:,2)]);

  
  rot_mat_rev=[cos(th) 0 sin(th); 0 1 0;-sin(th) 0 cos(th)]; % rot_mat_y
  XYZcr=rot_mat_rev*[xc,yc,0]';
  XYZc(sp,:)=bsxfun(@plus,XYZcr',scentroid);
  
  if DIAG
    figure(2), clf; hold on;
    plot3(XYZc(sp,1),XYZc(sp,2),XYZc(sp,3),'*r');
    plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'o');
    plot3(scentroid(1),scentroid(2),scentroid(3),'sg');
    grid on, axis equal;
    
    
    figure(1), clf, hold on;
    scatter3(XYZ0(:,1),XYZ0(:,2),XYZ0(:,3),'o','fill');
    scatter3(XYZr(:,1),XYZr(:,2),XYZr(:,3),'o','fill');
    axis equal
    grid on;
  end
  D(sp)=distance(XYZc(sp,:),scentroid);
end

save(['speaker_pos_' data_date '_circle_fit.mat'],'XYZc');

