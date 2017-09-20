% 2015 12 28  Plots for multi-freq composite clicks
% 2016 04 20  Separate this part out from multifreq_composite_click.m
% 2016 05 07  Plot for NIFTI poster
% 2016 05 07  Plot for ASA talk
% 2016 08 09  Plot for paper;
%             revised from multifreq_composite_click_cntr_20160808
% 2016 10 25  Update for version 1025 with out-of-bound points
% 2017 03 08  Re-processed all files due to errors in mic_bp beampattern
%             compensation and degF to degC conversion.

clear

usrn = getenv('username');
if strcmp(usrn,'Wu-Jung')
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
else
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
end

save_plot_opt = 1;

% Set up various paths
data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';

results_path = 'analysis_results_figs';
data_path = 'multifreq_composite_click_20170308_batall_bin10_th0';

ss = strsplit(data_path,'_');
threshold = str2double(ss{end}(3:end));  % only use averaged results if number of measurement > threshold
binsize = str2double(ss{end-1}(4:end));  % bin size in azimuth and elevation
bat_num = ss{end-2}(4:end);

data_file = sprintf('%s_merged_clicks.mat',data_path);
load(fullfile(data_base_path,results_path,data_path,data_file));

[~,script_name,~] = fileparts(mfilename('fullpath'));
results_path = 'analysis_results_figs';
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

cgrey = 200*ones(1,3)/255;
num_freq = length(param.freq_wanted);
num_freq_plot = num_freq-2;
colorset = jet(num_freq_plot);
contour_sm_len = 10;

fig_left = figure;
axesm eckert4
axis off
framem('fedgecolor',cgrey,'flonlimit',[-180 180]);
gridm('gcolor',cgrey,'glinestyle','-');
hold on

fig_right = figure;
axesm eckert4
axis off
framem('fedgecolor',cgrey,'flonlimit',[-180 180]);
gridm('gcolor',cgrey,'glinestyle','-');
hold on

save_fname_pre = sprintf('%s_bat%s_bin%d_th%d',script_name,bat_num,binsize,threshold);

for iF=2:num_freq-1
    color_iF = iF-1;
    
    xy_sm_left(:,1) = smooth(multifreq_3dB_contour.left{iF}(:,1),contour_sm_len);
    xy_sm_left(:,2) = smooth(multifreq_3dB_contour.left{iF}(:,2),contour_sm_len);
    xy_sm_left(isnan(xy_sm_left(:,1)),:) = NaN;
    
    xy_sm_right(:,1) = smooth(multifreq_3dB_contour.right{iF}(:,1),contour_sm_len);
    xy_sm_right(:,2) = smooth(multifreq_3dB_contour.right{iF}(:,2),contour_sm_len);
    xy_sm_right(isnan(xy_sm_right(:,1)),:) = NaN;
    
    figure(fig_left)  % left click
    plot(xy_sm_left(:,1),xy_sm_left(:,2),'linewidth',3,'color',colorset(color_iF,:));
    
    figure(fig_right)  % right click
    plot(xy_sm_right(:,1),xy_sm_right(:,2),'linewidth',3,'color',colorset(color_iF,:));
    
    clear xy_sm_*
end

figure(fig_left)
colormap(jet(num_freq_plot))
colorbar('Ticks',linspace(0+1/num_freq_plot/2,1-1/num_freq_plot/2,num_freq_plot),...
    'TickLabels',{num2str(param.freq_wanted(2:end-1)'/1e3)},'location','southoutside');
grid
tightmap
title(sprintf('Averaged left click, bat %s, th=%d, bin=%ddeg',...
    bat_num,binsize,threshold));

saveas(fig_left,fullfile(save_path,[save_fname_pre,'_left.fig']),'fig');
saveSameSize(fig_left,'file',fullfile(save_path,[save_fname_pre,'_left.png']),...
    'format','png','renderer','painters');


figure(fig_right)
colormap(jet(num_freq_plot))
colorbar('Ticks',linspace(0+1/num_freq_plot/2,1-1/num_freq_plot/2,num_freq_plot),...
    'TickLabels',{num2str(param.freq_wanted(2:end-1)'/1e3)},'location','southoutside');
grid
tightmap
title(sprintf('Averaged right click, bat %s, th=%d, bin=%ddeg',...
    bat_num,binsize,threshold));

saveas(fig_right,fullfile(save_path,[save_fname_pre,'_right.fig']),'fig');
saveSameSize(fig_right,'file',fullfile(save_path,[save_fname_pre,'_right.png']),...
    'format','png','renderer','painters');



