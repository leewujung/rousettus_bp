% 2015 12 28  Averaged elevation cut for individual clicks

clear

usrn = getenv('username');
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);

plot_opt = 1;
save_opt = 1;

% Load compiled rotated data
base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
compile_path = '20151224_all_click_rotate';
compile_file = 'all_click_rotated_34271.mat';
load(fullfile(base_path,compile_path,compile_file));

% Fix path with current username
if ~strcmp(usrn,'wlee76')
    if strfind(base_path, 'wlee76')
        base_path = regexprep(base_path,'wlee76',usrn);
    end
end

% Path/filename for saving data and figures
save_header = 'rousettus_34271';
save_path = fullfile(base_path,'20160209_mean_az_belt');  % set path for saving files
if ~exist(save_path,'dir')
    mkdir(save_path);
end


A.base_path = base_path;
A.compile_path = compile_path;
A.compile_file = compile_file;

% Frequency to be plotted
freq_wanted = 35e3;
A.freq_wanted  = freq_wanted;

% Elevation range to be averaged over
el_cut = 30;
el_range = -el_cut:3:el_cut;
az_range = -180:3:180;
[az_band,el_band] = meshgrid(az_range,el_range);
tot_area = range(el_range)*range(az_range);

% Get averaged elevation band
cut_az = cell(length(processed_files),1);
click_side = cell(length(processed_files),1);
for iB=1:length(processed_files)
    num_good_call = size(raw_meas.call_dB{iB},1);
    cut_az{iB} = nan(num_good_call,length(az_range));
    click_side{iB} = nan(num_good_call,1);
    for iC=1:num_good_call
        call_dB = raw_meas.call_dB{iB}(iC,:);
        call_dB_norm = call_dB-max(call_dB);
        az = shift_tilt_final.az{iB}(iC,:);
        el = shift_tilt_final.el{iB}(iC,:);
        ch_include_idx = raw_meas.ch_include_idx{iB}(iC,:) & ~isnan(az);
        vq_azel = rbfinterp([az_band(:)';el_band(:)'],...
                  rbfcreate([az(ch_include_idx);el(ch_include_idx)],call_dB_norm(ch_include_idx),...
                            'RBFFunction','multiquadrics'));
        vq_azel = reshape(vq_azel,size(az_band));
        
        % Set values outside of boundary to NaN
        azk = az(ch_include_idx);
        elk = el(ch_include_idx);
        k = boundary(azk',elk',0);  % outer boundary of all measured points
        [in,on] = inpolygon(az_band,el_band,azk(k),elk(k));
        in = in|on;
        vq_azel(~in) = NaN;
        
        % Get cut area
        az_band_cut = az_band(in);
        el_band_cut = el_band(in);
        kk = boundary(az_band_cut(:),el_band_cut(:),0);  % outer boundary of all measured points
        cut_area = polyarea(az_band_cut(kk),el_band_cut(kk));
        
%         if cut_area/tot_area>0.5  % criteria to include data
        cut_az{iB}(iC,:) = nanmean(vq_azel,1);
        click_side{iB}(iC) = raw_meas.click_side{iB}(iC);
%         end

    end
end

% Process left/right clicks separately
click_side_all = cell2mat(click_side);
cut_az_all = cell2mat(cut_az);

cut_az_right = cut_az_all(click_side_all==1,:);
cut_az_left = cut_az_all(click_side_all==0,:);

cut_az_right_mean = nanmean(cut_az_all(click_side_all==1,:));
cut_az_left_mean = nanmean(cut_az_all(click_side_all==0,:));

cut_az_right_std = nanstd(cut_az_all(click_side_all==1,:));
cut_az_left_std = nanstd(cut_az_all(click_side_all==0,:));

% Plot
tt = save_header;
tt(tt=='_') = ' ';

figure   % plot mean with all curves, separately
subplot(211)
plot(az_range,cut_az_all(click_side_all==1,:),'color',[1 1 1]*190/255)
hold on
hh = plot(az_range,cut_az_right_mean,'k','linewidth',1.5);
ylim([-40 5])
grid
ylabel('Relative intensity (dB)');
ll = legend(hh,'Right');
set(ll,'fontsize',12);
subplot(212)
plot(az_range,cut_az_all(click_side_all==0,:),'color',[1 1 1]*190/255)
hold on
hh = plot(az_range,cut_az_left_mean,'k','linewidth',1.5);
ylim([-40 5])
grid
ll = legend(hh,'Left');
set(ll,'fontsize',12);
xlabel('Azimuth (deg)');
ylabel('Relative intensity (dB)');
suptitle([tt,', el\_cut=',num2str(el_cut),'deg']);
if save_opt==1
    saveas(gcf,fullfile(save_path,[save_header,'_elcut_',num2str(el_cut),'deg_allcurves.png']),'png');
end


notnanidx_right = ~isnan(cut_az_right_mean) & ~isnan(cut_az_right_std);
notnanidx_left = ~isnan(cut_az_left_mean) & ~isnan(cut_az_left_std);


figure;  % Plot mean with std, separately
subplot(211)
patch([az_range(notnanidx_right),fliplr(az_range(notnanidx_right))],...
      [cut_az_right_mean(notnanidx_right)-cut_az_right_std(notnanidx_right),...
       fliplr(cut_az_right_mean(notnanidx_right)+cut_az_right_std(notnanidx_right))],...
       'r','edgecolor','none','facealpha',0.2);
hold on
hright = plot(az_range,cut_az_right_mean,'r','linewidth',2);
grid
ylim([-40 0])
ylabel('Relative intensity (dB)');
set(gca,'fontsize',12);
ll = legend(hright,'Right');
set(ll,'fontsize',12);
subplot(212)
patch([az_range(notnanidx_left),fliplr(az_range(notnanidx_left))],...
      [cut_az_left_mean(notnanidx_left)-cut_az_left_std(notnanidx_left),...
       fliplr(cut_az_left_mean(notnanidx_left)+cut_az_left_std(notnanidx_left))],...
       'b','edgecolor','none','facealpha',0.2);
hold on
hleft = plot(az_range,cut_az_left_mean,'b','linewidth',2);
grid
ylim([-40 0])
xlabel('Azimuth (deg)');
ylabel('Relative intensity (dB)');
set(gca,'fontsize',12);
ll = legend(hleft,'Left');
suptitle([tt,', el\_cut=',num2str(el_cut),'deg']);
if save_opt==1
    saveas(gcf,fullfile(save_path,[save_header,'_elcut_',num2str(el_cut),'deg_separate.png']),'png');
end


figure;  % Plot mean with std, overlapping
patch([az_range(notnanidx_right),fliplr(az_range(notnanidx_right))],...
      [cut_az_right_mean(notnanidx_right)-cut_az_right_std(notnanidx_right),...
       fliplr(cut_az_right_mean(notnanidx_right)+cut_az_right_std(notnanidx_right))],...
       'r','edgecolor','none','facealpha',0.2);
hold on
hright = plot(az_range,cut_az_right_mean,'r','linewidth',2);
patch([az_range(notnanidx_left),fliplr(az_range(notnanidx_left))],...
      [cut_az_left_mean(notnanidx_left)-cut_az_left_std(notnanidx_left),...
       fliplr(cut_az_left_mean(notnanidx_left)+cut_az_left_std(notnanidx_left))],...
       'b','edgecolor','none','facealpha',0.2);
hleft = plot(az_range,cut_az_left_mean,'b','linewidth',2);
grid
ll = legend([hright,hleft],'Right','Left');
set(ll,'fontsize',12);
xlabel('Azimuth (deg)');
ylabel('Relative intensity (dB)');
set(gca,'fontsize',12);
title([tt,', el\_cut=',num2str(el_cut),'deg']);
if save_opt==1
    saveas(gcf,fullfile(save_path,[save_header,'_elcut_',num2str(el_cut),'deg_combine.png']),'png');
end
