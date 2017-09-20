% 2015 12 28  Plots for multi-freq composite clicks
% 2016 04 20  Separate this part out from multifreq_composite_click.m
% 2016 05 07  Plot for NIFTI poster
% 2016 05 07  Plot for ASA talk
% 2016 08 08  Plot for paper
% 2016 08 13  Adapted to plot clustered avg bp

clear

usrn = getenv('username');
if strcmp(usrn,'Wu-Jung')
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
else
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
end

save_plot_opt = 1;

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
else
    data_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
    save_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
end

results_path = 'analysis_results_figs';
data_path = 'cluster_multifreq_composite_click_20160813_bat36134_bin10_th0';
save_path = fullfile(save_base_path,results_path,data_path);

ss = strsplit(data_path,'_');
threshold = str2double(ss{end}(3:end));  % only use averaged results if number of measurement > threshold
binsize = str2double(ss{end-1}(4:end));  % bin size in azimuth and elevation
bat_num = ss{end-2}(4:end);

knum = 3;
el_cut = 3;
freq_cluster = 35e3;

data_file = sprintf('%s_knum%d_elcut%02d_freqcls%02dkHz',...
    data_path,knum,el_cut,freq_cluster/1e3);

% Load clustering results
load(fullfile(data_base_path,results_path,data_path,data_file));

% Load elcut results
C = load(fullfile(data_base_path,results_path,param.elcut_path,param.elcut_file));


% Plot
cgrey = 200*ones(1,3)/255;
vq_norm_min = -27;
contour_vec = 0:-3:(floor(vq_norm_min/3)-1)*3;
cvec_min_idx = find(contour_vec-vq_norm_min<0,1,'first');
num_freq = length(param.freq_wanted);
iF = find(freq_cluster==param.freq_wanted);
% for iF=2:2:num_freq
fig = figure('units','normalized','outerposition',[0 0 1 1]);
%     fig = figure('position',[200,200,1200,450]);

idx_right = find(scatter_data.click_side_all==1);
idx_left = find(scatter_data.click_side_all==0);

az_max_right = scatter_data.az_max_all(scatter_data.click_side_all==1);
azq_max_right = scatter_data.azq_max_all(scatter_data.click_side_all==1);
el_max_right = scatter_data.el_max_all(scatter_data.click_side_all==1);
elq_max_right = scatter_data.elq_max_all(scatter_data.click_side_all==1);

az_max_left = scatter_data.az_max_all(scatter_data.click_side_all==0);
azq_max_left = scatter_data.azq_max_all(scatter_data.click_side_all==0);
el_max_left = scatter_data.el_max_all(scatter_data.click_side_all==0);
elq_max_left = scatter_data.elq_max_all(scatter_data.click_side_all==0);

for iK=1:knum
    % Mean curves
    subplot(3,2,1);
    plot(C.az_range,nanmean(C.cut_az_left(param.idx_cluster_left==iK,:),1),'linewidth',2);
    hold on
    subplot(3,2,2);
    plot(C.az_range,nanmean(C.cut_az_right(param.idx_cluster_right==iK,:),1),'linewidth',2);
    hold on
    
    idx_right_k = idx_right(param.idx_cluster_right==iK);
    idx_left_k = idx_left(param.idx_cluster_left==iK);
    
    subplot(3,2,3);
    plot(az_max_left(idx_left_k),el_max_left(idx_left_k),'.');
    hold on
    
    subplot(3,2,4);
    plot(az_max_right(idx_right_k),el_max_right(idx_right_k),'.');
    hold on
    
    subplot(3,2,5);
    plot(azq_max_left(idx_left_k),elq_max_left(idx_left_k),'.');
    hold on
    
    subplot(3,2,6);
    plot(azq_max_right(idx_right_k),elq_max_right(idx_right_k),'.');
    hold on
    
end % cluster loop

subplot(3,2,1);
xlabel('Azimuth (deg)'); ylabel('Norm''ed energy');
xlim([-180 180]); ylim([-30 3]); grid on
legend({num2str([1:knum]')});
subplot(3,2,2);
xlabel('Azimuth (deg)'); ylabel('Norm''ed energy');
xlim([-180 180]); ylim([-30 3]); grid on
legend({num2str([1:knum]')});

subplot(3,2,3);
xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
xlim([-180 180]); ylim([-90 90]); grid on
legend({num2str([1:knum]')});
title('Left, az/el raw');
subplot(3,2,4);
xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
xlim([-180 180]); ylim([-90 90]); grid on
legend({num2str([1:knum]')});
title('Right, az/el raw');

subplot(3,2,5);
xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
xlim([-180 180]); ylim([-90 90]); grid on
legend({num2str([1:knum]')});
title('Left, az/el interp');
subplot(3,2,6);
xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
xlim([-180 180]); ylim([-90 90]); grid on
legend({num2str([1:knum]')});
title('Right, az/el interp');

suptitle(sprintf('Az/el, clustered at %dkHz, bat %s, %dkHz th=%d, bin=%d, knum=%d',...
    param.freq_used_for_cluster/1e3,bat_num,param.freq_wanted(iF)/1e3,...
    param.composite_threshold,param.composite_binsize,knum));

% Save figures
save_fname = sprintf('%s_knum%d_%dkHz_azel_distr',...
    data_file,knum,param.freq_wanted(iF)/1e3);
%     saveas(fig,fullfile(save_path,[save_fname,'.fig']),'fig');
saveSameSize(fig,'file',fullfile(save_path,[save_fname,'.png']),...
    'format','png','renderer','painters');

close(fig)

% end % freq loop


