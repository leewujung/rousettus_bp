% 2015 12 28  Plots for multi-freq composite clicks
% 2016 04 20  Separate this part out from multifreq_composite_click.m
% 2016 05 07  Plot for NIFTI poster
% 2016 05 07  Plot for ASA talk
% 2016 08 08  Plot for paper

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
data_path = 'multifreq_composite_click_20160808_bat34271_bin10_th0';
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
    
    fig = figure('position',[200,200,1200,450]);
    suptitle(sprintf('Averaged right clicks, bat %s, %d kHz, th=%d, bin=%ddeg',...
        bat_num,param.freq_wanted(iF)/1e3,param.composite_threshold,param.composite_binsize));

    subplot(121)    % left click
    axesm eckert4
    framem('fedgecolor',cgrey,'flonlimit',[-180 180]);
    gridm('gcolor',cgrey,'glinestyle','-');
    axis off
    contourfm(averaged_composite.left.interp(iF).elq_avg,...
        averaged_composite.left.interp(iF).azq_avg,...
        averaged_composite.left.interp(iF).vq_norm_avg,...
        contour_vec(2:cvec_min_idx),...
        'fill','on','linecolor','w');  % don't plot 0 dB contour
    colorbar('southoutside','ticks',fliplr(contour_vec(1:cvec_min_idx)));
    colormap(parula(cvec_min_idx-1));
    caxis([contour_vec(cvec_min_idx) 0]);
    tightmap
    title('Left');

    subplot(122);  % right click
    axesm eckert4
    framem('fedgecolor',cgrey,'flonlimit',[-180 180]);
    gridm('gcolor',cgrey,'glinestyle','-');
    axis off
    contourfm(averaged_composite.right.interp(iF).elq_avg,...
        averaged_composite.right.interp(iF).azq_avg,...
        averaged_composite.right.interp(iF).vq_norm_avg,...
        contour_vec(2:cvec_min_idx),...
        'fill','on','linecolor','w');  % don't plot 0 dB contour
    colorbar('southoutside','ticks',fliplr(contour_vec(1:cvec_min_idx)));
    colormap(parula(cvec_min_idx-1));
    caxis([contour_vec(cvec_min_idx) 0]);
    tightmap
    title('Right');
    
    % Save figures
    save_fname = sprintf('%s_%02dkHz_avg_bp',data_path,param.freq_wanted(iF)/1e3);
    saveas(fig,fullfile(save_path,[save_fname,'.fig']),'fig');
    saveSameSize(fig,'file',fullfile(save_path,[save_fname,'.png']),...
        'format','png','renderer','painters');
    
    close(fig)
end



