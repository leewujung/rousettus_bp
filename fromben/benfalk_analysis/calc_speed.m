function speed=calc_speed(bat_wo_wb,fvideo,DIAG)



speed=sqrt(sum( diff( bat_wo_wb,[], 1 ).^2 , 2 ) ) * fvideo ;



% val_ind = find(~isnan(centroid(:,1)));

% interp_xval = smooth(interp1(val_ind,centroid(val_ind,1),1:length(centroid)),5);
% interp_yval = smooth(interp1(val_ind,centroid(val_ind,2),1:length(centroid)),5);
% interp_zval = smooth(interp1(val_ind,centroid(val_ind,3),1:length(centroid)),5);
% 
% interp_xval(interp_xval==0)=nan;
% interp_yval(interp_yval==0)=nan;
% interp_zval(interp_zval==0)=nan;
% 
% sm_centroid = [interp_xval interp_yval interp_zval];
% 
% speed=sqrt(sum( diff( sm_centroid ).^2 , 2 ) ) * fvideo ;
% speed(speed==0)=nan;
% 
% nan_indx=find(isnan(centroid(:,1)));
% sm_centroid(nan_indx,:)=nan;
% speed(nan_indx,:)=nan;

if nargin >= 3 && DIAG  
  figure(14); clf;
  plot(speed);
  title('speed (m/s)')
end