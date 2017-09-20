function [cam_p,cam_o,cam_n]=extract_cam_positions(file)

s=xml2struct(file);

cam_p=nan(length(s.Cameras.Camera),3);
cam_o=nan(length(s.Cameras.Camera),4);
cam_n=nan(length(s.Cameras.Camera),1);
for k=1:length(s.Cameras.Camera)
  C=s.Cameras.Camera{k};
  cam_p(k,:)=sscanf(C.KeyFrames.KeyFrame.Attributes.POSITION,'%f %f %f');
  
  %orientation is in quaternions with the w element last
  %to extract orientation direction, convert to Rot Mat and use 3rd row
  %R = quat2rot(q([4 1:3]));
  %E=R(3,:);
  cam_o(k,:)=sscanf(C.KeyFrames.KeyFrame.Attributes.ORIENTATION,'%f %f %f %f');
  
  cam_n(k)=sscanf(C.Attributes.USERID,'%f')+1;
end