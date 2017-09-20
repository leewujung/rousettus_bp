%match mics with rotation and translation
clear

thresh_dist=.01;

mic_dir='..\mic_pos\';

date_data='20150911';
date_data2='20150910';

trl1='_S1_01'; %one you want to fix
trl2='_S1_04'; %one you want to pull from

mics1=load([mic_dir 'mic_pos_' date_data trl1 '.mat']);
mics2=load([mic_dir 'mic_pos_' date_data2 trl2 '.mat']);

order1=load([mic_dir 'mic_order_' date_data trl1 '.mat']);
order2=load([mic_dir 'mic_order_' date_data2 trl2 '.mat']);



% if not using calib mics...
if datenum(date_data,'yyyymmdd') < 736200 %20150824
  no_calib=1;
else % using calib mics:
  no_calib=0;
end

if no_calib
  tip_indx=(1:3:32*3)+2;
  base_indx=[ (tip_indx(1:32)-1)  ;...
    (tip_indx(1:32)-2)] ;
else
  tip_indx=[(1:3:32*3)+2 32*3+2 32*3+4];
  base_indx=[ (tip_indx(1:32)-1) tip_indx([33 34]) ;...
    (tip_indx(1:32)-2) tip_indx([33 34])-1] ;
end

tips1_unorder=mics1.mic_pos(tip_indx,:);
tips2_unorder=mics2.mic_pos(tip_indx,:);
base1_unorder=mics1.mic_pos(base_indx(1,:),:);
base2_unorder=mics2.mic_pos(base_indx(1,:),:);
mid1_unorder=mics1.mic_pos(base_indx(2,:),:);
mid2_unorder=mics2.mic_pos(base_indx(2,:),:);

order1=order1.mic_order;
order2=order2.mic_order;

%ordering
[tips1,tips2,base1,base2,mid1,mid2]=deal(nan(size(tips1_unorder)));
tips1(isfinite(order1),:)=tips1_unorder(order1(isfinite(order1)),:);
tips2(isfinite(order2),:)=tips2_unorder(order2(isfinite(order2)),:);
base1(isfinite(order1),:)=base1_unorder(order1(isfinite(order1)),:);
base2(isfinite(order2),:)=base2_unorder(order2(isfinite(order2)),:);
mid1(isfinite(order1),:)=mid1_unorder(order1(isfinite(order1)),:);
mid2(isfinite(order2),:)=mid2_unorder(order2(isfinite(order2)),:);


%using tips
% indx=~isnan(tips1(:,1)) & ~isnan(tips2(:,1));
% A=tips1(indx,:);
% B=tips2(indx,:);
% tips2trans=(R*tips2'+repmat(T,1,size(tips2,1)))';

%using base
indx=~isnan(base1(:,1)) & ~isnan(base2(:,1)) & ~isnan(mid1(:,1))...
   & ~isnan(mid2(:,1));
A=[base1(indx,:); mid1(indx,:)];
B=[base2(indx,:); mid2(indx,:)];

[R,T]=rigid_transform_3D(B,A);

base2trans=(R*base2'+repmat(T,1,size(base2,1)))';
mid2trans=(R*mid2'+repmat(T,1,size(mid2,1)))';
tips2trans=(R*tips2'+repmat(T,1,size(tips2,1)))';

MDbase=distance(base1,base2trans);
MDmid=distance(mid1,mid2trans);

%indx with low residuals
indx = MDbase < thresh_dist & MDmid < thresh_dist;
if sum(indx) > 10
  A=[base1(indx,:); mid1(indx,:)];
  B=[base2(indx,:); mid2(indx,:)];

  [R,T]=rigid_transform_3D(B,A);

  base2trans=(R*base2'+repmat(T,1,size(base2,1)))';
  mid2trans=(R*mid2'+repmat(T,1,size(mid2,1)))';
  tips2trans=(R*tips2'+repmat(T,1,size(tips2,1)))';

  MDbase=distance(base1,base2trans);
  MDmid=distance(mid1,mid2trans);
end

fprintf('\nBase difference after rotating and translating\n')
fprintf('Mean: %2.4f\n',nanmean(MDbase));
fprintf('std: %2.4f\n',nanstd(MDbase));
fprintf('Max: %2.4f\n',max(MDbase));
fprintf('Min: %2.4f\n',min(MDbase));

fprintf('\nMid difference after rotating and translating\n')
fprintf('Mean: %2.4f\n',nanmean(MDmid));
fprintf('std: %2.4f\n',nanstd(MDmid));
fprintf('Max: %2.4f\n',max(MDmid));
fprintf('Min: %2.4f\n',min(MDmid));

missingbase1=isnan(base1(:,1));
missingmid1=isnan(mid1(:,1));
missingtip1=isnan(tips1(:,1));

base1rep=base1;
mid1rep=mid1;
tips1rep=tips1;
base1rep(missingbase1,:)=base2trans(missingbase1,:);
mid1rep(missingmid1,:)=mid2trans(missingmid1,:);
tips1rep(missingtip1,:)=tips2trans(missingtip1,:);

tips=tips1rep;
mic_vec=base1rep-mid1rep;

save(['..\mic_pos\mics_ch_order_' date_data trl1(1:6) '_merged.mat'],...
  'tips','mic_vec');

figure(3), clf; set(gcf,'pos',[10 300 912 700]), hold on;
scatter3(tips1(:,3),tips1(:,1),tips1(:,2),'k')
text(tips1(:,3),tips1(:,1),tips1(:,2),num2str((1:length(tips1(:,1)))'))

hold on

scatter3(tips2trans(:,3),tips2trans(:,1),tips2trans(:,2),'+')
text(tips2trans(:,3),tips2trans(:,1),tips2trans(:,2),num2str((1:length(tips2trans(:,1)))'),...
  'color','r')

axis equal, grid on  
view(-69.2,22)