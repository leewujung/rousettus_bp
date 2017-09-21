% 2015 12 28  Plots for multi-freq composite clicks
% 2016 04 20  Separate this part out from multifreq_composite_click.m
% 2016 05 07  Plot for NIFTI poster
% 2016 05 07  Plot for ASA talk
% 2016 08 08  Plot for paper
% 2017 03 08  Re-processed all files due to errors in mic_bp beampattern
%             compensation and degF to degC conversion.
% 2017 09 20  -- Plot beam center on the composite beampattern
%             -- Check location of the center of best-fitting ellipse

clear
if isunix
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/rbfinterp_v1.2');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/EllipseDirectFit');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/saveSameSize');

    data_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
    save_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
else
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    
    data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
end


save_plot_opt = 1;

% Set up various paths
results_path = 'analysis_results_figs';
data_path = 'multifreq_composite_click_20170308_batall_bin10_th0';

ss = strsplit(data_path,'_');
threshold = str2double(ss{end}(3:end));  % only use averaged results if number of measurement > threshold
binsize = str2double(ss{end-1}(4:end));  % bin size in azimuth and elevation
bat_num = ss{end-2}(4:end);

data_file = sprintf('%s_merged_clicks.mat',data_path);
load(fullfile(data_base_path,results_path,data_path,data_file));

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end
save_fname_pre = sprintf('%s_bat%s_bin%d_th%d',...
    script_name,bat_num,binsize,threshold);

cgrey = 200*ones(1,3)/255;
num_freq = length(param.freq_wanted);
for iF=2:2:num_freq
    
    vq_norm_min = -27;
    contour_vec = 0:-3:(floor(vq_norm_min/3)-1)*3;
    %contour_vec = 0:-1:(floor(vq_norm_min/3)-1)*3;
    cvec_min_idx = find(contour_vec-vq_norm_min<0,1,'first');

    % Find max beam energy point
    % --- left composite click
    xx = averaged_composite.left.interp(iF).vq_norm_avg(:);
    [left_max_val,left_max_idx] = max(xx);
    left_max_el = averaged_composite.left.interp(iF).elq_avg(left_max_idx);
    left_max_az = averaged_composite.left.interp(iF).azq_avg(left_max_idx);
    xx(isnan(xx)) = -inf;
    [~,left_sort_idx] = sort(xx,'descend');
    ii = xx(left_sort_idx)>-1;
    left_top_el = mean(averaged_composite.left.interp(iF).elq_avg(left_sort_idx(ii)));
    left_top_az = mean(averaged_composite.left.interp(iF).azq_avg(left_sort_idx(ii)));
    % --- right composite click
    xx = averaged_composite.right.interp(iF).vq_norm_avg(:);
    [right_max_val,right_max_idx] = max(xx);
    right_max_el = averaged_composite.right.interp(iF).elq_avg(right_max_idx);
    right_max_az = averaged_composite.right.interp(iF).azq_avg(right_max_idx);
    xx(isnan(xx)) = -inf;
    [~,right_sort_idx] = sort(xx,'descend');
    ii = xx(right_sort_idx)>-1;
    right_top_el = mean(averaged_composite.right.interp(iF).elq_avg(right_sort_idx(ii)));
    right_top_az = mean(averaged_composite.right.interp(iF).azq_avg(right_sort_idx(ii)));

    % Fit ellipse
    map_proj = 'eckert4';
    mstruct = defaultm(map_proj);
    mstruct = defaultm(mstruct);
    % --- left composite click
    [left_raw,left_rot_max,left_rot_elpctr,left_rot_elpctr_tilt] = ...
        shift_rotate_bp_composite(averaged_composite.left.bin(iF).az_avg,...
                                  averaged_composite.left.bin(iF).el_avg,...
                                  averaged_composite.left.bin(iF).avg_call_dB,...
                                  map_proj,0.005);
    [left_el_ectr,left_az_ectr] = minvtran(mstruct,left_rot_max.E.x0,left_rot_max.E.y0);  % inverse map projection
    [left_el_ectr_r,left_az_ectr_r] = rotatem(left_el_ectr,left_az_ectr,...
                                    [left_max_el,left_max_az],...
                                    'inverse','degrees');
    % --- right composite click
    [right_raw,right_rot_max,right_rot_elpctr,right_rot_elpctr_tilt] = ...
        shift_rotate_bp_composite(averaged_composite.right.bin(iF).az_avg,...
                                  averaged_composite.right.bin(iF).el_avg,...
                                  averaged_composite.right.bin(iF).avg_call_dB,...
                                  map_proj,0.005);
    [right_el_ectr,right_az_ectr] = minvtran(mstruct,right_rot_max.E.x0,right_rot_max.E.y0);  % inverse map projection
    [right_el_ectr_r,right_az_ectr_r] = rotatem(right_el_ectr,right_az_ectr,...
                                    [right_max_el,right_max_az],...
                                    'inverse','degrees');

    
    % Plot left click
    figure
    axesm eckert4
    gridm('gcolor',cgrey,'glinestyle','-');
    framem('fedgecolor',cgrey,'flonlimit',[-180 180]);
    axis off
    contourfm(averaged_composite.left.interp(iF).elq_avg,...
        averaged_composite.left.interp(iF).azq_avg,...
        averaged_composite.left.interp(iF).vq_norm_avg,...
        contour_vec(2:cvec_min_idx),...
        'fill','on','linecolor','w');  % don't plot 0 dB contour
    
    % Plot max beam energy point
    plotm(left_max_el,left_max_az,'rx','markersize',8,'linewidth',2)

    % Plot dB>-1 beam energy point
    plotm(left_top_el,left_top_az,'r^','markersize',8,'linewidth',2)

    % location of center of best-fitting ellipse
    plotm(left_el_ectr_r,left_az_ectr_r,'ro','markersize',8,'linewidth',2)

    colorbar('southoutside','ticks',fliplr(contour_vec(1:cvec_min_idx)));
    colormap(parula(cvec_min_idx-1));
    caxis([contour_vec(cvec_min_idx) 0]);
    tightmap
    title(sprintf('Averaged left clicks, bat %s, %d kHz, th=%d, bin=%ddeg',...
        bat_num,param.freq_wanted(iF)/1e3,param.composite_threshold,param.composite_binsize));

    % Save figures
    save_fname = sprintf('%s_%02dkHz_avg_bp_left',...
        save_fname_pre,param.freq_wanted(iF)/1e3);
    saveas(gcf,fullfile(save_path,[save_fname,'.fig']),'fig');
    saveSameSize(gcf,'file',fullfile(save_path,[save_fname,'.png']),...
        'format','png','renderer','painters');
    close
    
    % Plot right click
    figure   
    axesm eckert4
    gridm('gcolor',cgrey,'glinestyle','-');
    framem('fedgecolor',cgrey,'flonlimit',[-180 180]);
    axis off
    contourfm(averaged_composite.right.interp(iF).elq_avg,...
        averaged_composite.right.interp(iF).azq_avg,...
        averaged_composite.right.interp(iF).vq_norm_avg,...
        contour_vec(2:cvec_min_idx),...
        'fill','on','linecolor','w');  % don't plot 0 dB contour
    % Plot max beam energy point
    plotm(right_max_el,right_max_az,'rx','markersize',8,'linewidth',2)

    % Plot dB>-1 beam energy point
    plotm(right_top_el,right_top_az,'r^','markersize',8,'linewidth',2)

    % location of center of best-fitting ellipse
    plotm(right_el_ectr_r,right_az_ectr_r,'ro','markersize',8,'linewidth',2)

    colorbar('southoutside','ticks',fliplr(contour_vec(1:cvec_min_idx)));
    colormap(parula(cvec_min_idx-1));
    caxis([contour_vec(cvec_min_idx) 0]);
    tightmap
    title(sprintf('Averaged right clicks, bat %s, %d kHz, th=%d, bin=%ddeg',...
        bat_num,param.freq_wanted(iF)/1e3,param.composite_threshold,param.composite_binsize));
    
    % Save figures
    save_fname = sprintf('%s_%02dkHz_avg_bp_right',...
        save_fname_pre,param.freq_wanted(iF)/1e3);
    saveas(gcf,fullfile(save_path,[save_fname,'.fig']),'fig');
    saveSameSize(gcf,'file',fullfile(save_path,[save_fname,'.png']),...
        'format','png','renderer','painters');
    close
end



