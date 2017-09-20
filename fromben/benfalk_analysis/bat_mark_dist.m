%bat marker distances
clear;

bat_pos_dir='..\bat_pos\';

files=dir([bat_pos_dir '*.mat']);
fnames={files.name};

%load in all files within a bat and within a day
%unique species
%unique days
%unique bands
FF=cellfun(@(c) strsplit(c,'_'),fnames,'uniformoutput',0);
all_sp=cellfun(@(c) c{1},FF,'uniformoutput',0);
un_sp=unique(all_sp);

all_day=cellfun(@(c) c{2},FF,'uniformoutput',0);
un_day=unique(all_day);

all_band=cellfun(@(c) c{3},FF,'uniformoutput',0);
un_band=unique(all_band);


for ss=1:length(un_sp)
  indx_sp=strcmp(all_sp,un_sp{ss});
  for dd=1:length(un_day)
    indx_day=strcmp(all_day,un_day{dd});
    for bb=1:length(un_band)
      indx_band=strcmp(all_band,un_band{bb});
      
      itrials = indx_sp & indx_day & indx_band;
      %calculate the mean distance for each of the markers for every frame for
      %all the trials
      
      if isempty(find(itrials, 1))
        continue
      end
      
      [DD12,DD13,DD23]=deal([]);
      for k=find(itrials)
        load([bat_pos_dir fnames{k}]);
        
        frames_w_alldata=find(isfinite(bat_pos{1}(:,1)) & isfinite(bat_pos{2}(:,1))...
          & isfinite(bat_pos{3}(:,1)));
        
        DD12=[DD12; distance(bat_pos{1}(frames_w_alldata,:),bat_pos{2}(frames_w_alldata,:))];
        DD13=[DD13; distance(bat_pos{1}(frames_w_alldata,:),bat_pos{3}(frames_w_alldata,:))];
        DD23=[DD23; distance(bat_pos{2}(frames_w_alldata,:),bat_pos{3}(frames_w_alldata,:))];
      end
      
      
      
      %get a mean and exclude outliers from that mean
      figure(1); clf; histogram(DD12)
      hold on;
      thresh = 2*std(DD12);
      a=axis;
      plot([-thresh+mean(DD12) -thresh+mean(DD12)],[0 a(4)]);
      plot([thresh+mean(DD12) thresh+mean(DD12)],[0 a(4)]);
      
      figure(2); clf; histogram(DD13)
      hold on;
      thresh = 2*std(DD13);
      a=axis;
      plot([-thresh+mean(DD13) -thresh+mean(DD13)],[0 a(4)]);
      plot([thresh+mean(DD13) thresh+mean(DD13)],[0 a(4)]);
      
      figure(3); clf; histogram(DD23)
      hold on;
      thresh = 2*std(DD23);
      a=axis;
      plot([-thresh+mean(DD23) -thresh+mean(DD23)],[0 a(4)]);
      plot([thresh+mean(DD23) thresh+mean(DD23)],[0 a(4)]);
      
      %go back to each trial and exclude the points whose distances are too far
      %above the mean distance for that day and bat
    end
  end
end



[MD12,MD13,MD23]=deal(nan(size(files)));
for ff=1:length(files);
  load([bat_pos_dir files(ff).name])

  frames_w_alldata=find(isfinite(bat_pos{1}(:,1)) & isfinite(bat_pos{2}(:,1))...
    & isfinite(bat_pos{3}(:,1)));

  DD12=distance(bat_pos{1}(frames_w_alldata,:),bat_pos{2}(frames_w_alldata,:));
  DD13=distance(bat_pos{1}(frames_w_alldata,:),bat_pos{3}(frames_w_alldata,:));
  DD23=distance(bat_pos{2}(frames_w_alldata,:),bat_pos{3}(frames_w_alldata,:));

  MD12(ff)=mean(DD12);
  MD13(ff)=mean(DD13);
  MD23(ff)=mean(DD23);
  
%   figure(1), clf;
%   histogram(DD12,mean(DD12) - .005 : .001 : mean(DD12)+.005);
% 
%   figure(2), clf;
%   histogram(DD13,mean(DD13) - .005 : .001 : mean(DD13)+.005);
end
mean(MD12)
mean(MD13)
mean(MD23)