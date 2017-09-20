%match_mics
%load in known ordered position
%match to unknown position based on nearest neighbor

clear;

date_data='20150824';
trl_num='01';

if datenum(date_data,'yyyymmdd') < 736204 %20150828 %
  %2D
  ordered_mic_data=load('..\mic_pos\mic_pos_20150824_S1_03.mat');
  %load in the order of the ordered data with the channel mappings on array
  order_of_ordered_mics=load('..\mic_pos\mic_order_20150824_S1_03.mat');
else
  %cross
  ordered_mic_data=load(['..\mic_pos\mic_pos_20150910_S1_04.mat']);
  %load in the order of the ordered data with the channel mappings on array
  order_of_ordered_mics=load('..\mic_pos\mic_order_20150910_S1_04.mat');
end

% if not using calib mics...
if datenum(date_data,'yyyymmdd') < 736200 %20150824
  no_calib=1;
else % using calib mics:
  no_calib=0;
end

load(['..\mic_pos\mic_pos_' date_data '_S1_' trl_num '.mat'])

if no_calib
  tip_indx=(1:3:32*3)+2;
else
  tip_indx=[(1:3:32*3)+2 32*3+2 32*3+4];
end
tips=mic_pos(tip_indx,:);
tips_ord=ordered_mic_data.mic_pos(tip_indx,:);

for gg=1:3 %running it a few times should get the right offsets eventually
  offset=[];
  for k=1:3
    x=[];
    x(:,1)=tips(:,k);
    x(:,2)=tips_ord(:,k);
    offset(k)=fminsearch(@(a) dist_minimizer(a,x),0);
    tips_ord(:,k)=tips_ord(:,k)-offset(k);
  end
end

ii=nan(length(tips(:,1)),1);
for mm=1:length(tips)
  if isfinite(tips(mm,1))
    DD=distance(tips(mm,:),tips_ord);
    [~,ii(mm)]=min(DD);
  end
end

figure(2); clf,
set(gcf,'pos',[10 100 612 900]), clf; hold on;
scatter3(tips(:,3),tips(:,1),tips(:,2))
text(tips(:,3),tips(:,1),tips(:,2),...
  num2str((1:length(tips(:,1)))'),'color','r')
scatter3(tips_ord(:,3),tips_ord(:,1),tips_ord(:,2))
% text(tips_ord(:,3),tips_ord(:,1),tips_ord(:,2),...
%   num2str((1:length(tips_ord(:,1)))'))
axis equal, grid on;
view(-57,20)


%save order using new index
[~,mic_order]=ismember(order_of_ordered_mics.mic_order,ii);
mic_order(mic_order==0)=nan;
if length(unique(mic_order))==length(mic_order)
  if exist(['..\mic_pos\mic_order_' date_data '_S1_' trl_num '.mat'],'file')
    MO=load(['..\mic_pos\mic_order_' date_data '_S1_' trl_num '.mat']);
    mic_order1=MO.mic_order;
  end
  if ~exist('mic_order1','var') || ~isequaln(mic_order,mic_order1)
    save(['..\mic_pos\mic_order_' date_data '_S1_' trl_num '.mat'],'mic_order');
  end
else
  disp('there are duplicates')
end