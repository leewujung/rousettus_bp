% mics_2_ch_order
%need the mic order from assign_mic_order.m and the mic positions
%extracted from extract_mics

clear;

date_data='20150824';
trl_num='01';


load(['..\mic_pos\mic_order_' date_data '_S1_' trl_num '.mat']) %mic_order, assign_mic_order.m
load(['..\mic_pos\mic_pos_' date_data '_S1_' trl_num '.mat']) %mic_pos, extract_mics.m


% if not using calib mics...
if datenum(date_data,'yyyymmdd') < 736200 %20150824
  no_calib=1;
else % using calib mics:
  no_calib=0;
end

if no_calib
  tip_indx=(1:3:32*3)+2;
  vec_indx=[ (tip_indx(1:32)-1)  ;...
    (tip_indx(1:32)-2)] ;
else
  tip_indx=[(1:3:32*3)+2 32*3+2 32*3+4];
  vec_indx=[ (tip_indx(1:32)-1) tip_indx([33 34]) ;...
    (tip_indx(1:32)-2) tip_indx([33 34])-1] ;
end

tips1=mic_pos(tip_indx,:);
mic_vec1=mic_pos(vec_indx(1,:),:) - mic_pos(vec_indx(2,:),:);

tips=nan(size(tips1));
mic_vec=nan(size(tips1));
tips(isfinite(mic_order),:)=tips1(mic_order(isfinite(mic_order)),:);
mic_vec(isfinite(mic_order),:)=mic_vec1(mic_order(isfinite(mic_order)),:);

% figure(1), clf,
% plot3([zeros(length(mic_vec),1)'; mic_vec(:,3)'],...
%   [zeros(length(mic_vec),1)'; mic_vec(:,1)'],...
%   [zeros(length(mic_vec),1)'; mic_vec(:,2)'],'k')
% grid on; axis equal
% title('mic vectors')
% 
% figure(2), clf,
% plot3(tips(:,3),tips(:,1),tips(:,2),'ok')
% hold on;
% plot3(0,0,0,'or')
% plot3([zeros(length(mic_vec),1)'; tips(:,3)'],...
%   [zeros(length(mic_vec),1)'; tips(:,1)'],...
%   [zeros(length(mic_vec),1)'; tips(:,2)'],'k')
% grid on; axis equal
% title('Tips')


figure(3), clf; set(gcf,'pos',[10 300 912 700]), hold on;
scatter3(tips(:,3),tips(:,1),tips(:,2),'.k')
text(tips(:,3),tips(:,1),tips(:,2),num2str((1:length(tips(:,1)))'),'color','r')
view(-69.2,22); axis equal, grid on  



save(['..\mic_pos\mics_ch_order_' date_data '_S1_' trl_num '.mat'],...
  'tips','mic_vec');






