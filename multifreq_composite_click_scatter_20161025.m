% 2015 12 28  Plots for multi-freq composite clicks
% 2016 04 20  Separate this part out from multifreq_composite_click.m
% 2016 05 07  Plot for NIFTI poster
% 2016 05 07  Plot for ASA talk
% 2016 08 08  Plot for paper
% 2016 10 25  Update for version 1025 with out-of-bound points

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
    data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
end

results_path = 'analysis_results_figs';
data_path = 'multifreq_composite_click_20161025_batall_bin10_th0';
save_path = fullfile(save_base_path,results_path,data_path);

ss = strsplit(data_path,'_');
threshold = str2double(ss{end}(3:end));  % only use averaged results if number of measurement > threshold
binsize = str2double(ss{end-1}(4:end));  % bin size in azimuth and elevation
bat_num = ss{end-2}(4:end);

data_file = sprintf('%s_merged_clicks.mat',data_path);

load(fullfile(data_base_path,results_path,data_path,data_file));

cgrey = 200*ones(1,3)/255;
num_freq = length(param.freq_wanted);
for iF=1:num_freq
    
    vq_norm_min = -27;
    contour_vec = 0:-3:(floor(vq_norm_min/3)-1)*3;
    cvec_min_idx = find(contour_vec-vq_norm_min<0,1,'first');

    % all scatter data points
    call_dB_merge = reshape(squeeze(scatter_data.call_dB_all(iF,:,:)),[],1);
    click_side_merge = reshape(scatter_data.click_side_all,[],1);
    az_elpctr_tilt_merge = reshape(scatter_data.az_elpctr_tilt_all,[],1);
    el_elpctr_tilt_merge = reshape(scatter_data.el_elpctr_tilt_all,[],1);
    
    % left/right click
    idx_right_good = find(click_side_merge==1 & ~isnan(call_dB_merge) &...
                          ~isnan(az_elpctr_tilt_merge));
    call_dB_right = call_dB_merge(idx_right_good);
    az_right = az_elpctr_tilt_merge(idx_right_good);
    el_right = el_elpctr_tilt_merge(idx_right_good);
    
    idx_left_good = find(click_side_merge==0 & ~isnan(call_dB_merge) &...
                         ~isnan(az_elpctr_tilt_merge));
    call_dB_left = call_dB_merge(idx_left_good);
    az_left = az_elpctr_tilt_merge(idx_left_good);
    el_left = el_elpctr_tilt_merge(idx_left_good);
    
    fig = figure('position',[200,200,1200,450]);
    suptitle(sprintf('Averaged right clicks, bat %s, %d kHz, th=%d, bin=%ddeg',...
        bat_num,param.freq_wanted(iF)/1e3,param.composite_threshold,param.composite_binsize));

    subplot(121)    % left click
    axesm eckert4
    framem('fedgecolor',cgrey,'flonlimit',[-180 180]);
    gridm('gcolor',cgrey,'glinestyle','-');
    axis off
    scatterm(el_left,az_left,10,call_dB_left,'fill');
    colorbar('southoutside','ticks',fliplr(contour_vec(1:cvec_min_idx)));
%     colormap(parula(cvec_min_idx-1));
    caxis([contour_vec(cvec_min_idx) 0]);
    tightmap
    title('Left');

    subplot(122);  % right click
    axesm eckert4
    framem('fedgecolor',cgrey,'flonlimit',[-180 180]);
    gridm('gcolor',cgrey,'glinestyle','-');
    axis off
    scatterm(el_right,az_right,10,call_dB_right,'fill');
    colorbar('southoutside','ticks',fliplr(contour_vec(1:cvec_min_idx)));
%     colormap(parula(cvec_min_idx-1));
    caxis([contour_vec(cvec_min_idx) 0]);
    tightmap
    title('Right');
    
    % Save figures
    save_fname = sprintf('%s_%02dkHz_scatter',data_path,param.freq_wanted(iF)/1e3);
    saveas(fig,fullfile(save_path,[save_fname,'.fig']),'fig');
    saveSameSize(fig,'file',fullfile(save_path,[save_fname,'.png']),...
        'format','png','renderer','painters');

    close(fig)
end



