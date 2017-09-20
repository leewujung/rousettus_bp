% 2015 12 28  Plots for multi-freq composite clicks
% 2016 04 20  Separate this part out from multifreq_composite_click.m
% 2016 05 07  Plot for NIFTI poster
% 2016 05 07  Plot for ASA talk
% 2016 08 08  Plot for paper
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
if strcmp(usrn,'Wu-Jung')
    data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
else
    data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
end

results_path = 'analysis_results_figs';
data_path = 'multifreq_composite_click_20170308_batall_bin10_th3';
save_path = fullfile(save_base_path,results_path,data_path);

ss = strsplit(data_path,'_');
threshold = str2double(ss{end}(3:end));  % only use averaged results if number of measurement > threshold
binsize = str2double(ss{end-1}(4:end));  % bin size in azimuth and elevation
bat_num = ss{end-2}(4:end);

data_file = sprintf('%s_merged_clicks.mat',data_path);

load(fullfile(data_base_path,results_path,data_path,data_file));


cgrey = 200*ones(1,3)/255;
num_freq = length(param.freq_wanted);
colorset = jet(num_freq);
contour_sm_len = 10;

fig = figure('position',[200,200,1200,450]);
suptitle(sprintf('Averaged clicks, bat %s, th=%d, bin=%ddeg',...
    bat_num,binsize,threshold));

for iF=1:num_freq
    
    xy_sm_left(:,1) = smooth(multifreq_3dB_contour.left{iF}(:,1),contour_sm_len);
    xy_sm_left(:,2) = smooth(multifreq_3dB_contour.left{iF}(:,2),contour_sm_len);
    xy_sm_left(isnan(xy_sm_left(:,1)),:) = NaN;
    
    xy_sm_right(:,1) = smooth(multifreq_3dB_contour.right{iF}(:,1),contour_sm_len);
    xy_sm_right(:,2) = smooth(multifreq_3dB_contour.right{iF}(:,2),contour_sm_len);
    xy_sm_right(isnan(xy_sm_right(:,1)),:) = NaN;
    
    figure(fig)
    subplot(121)  % left click
    if iF==1
        axesm eckert4
        axis off
        framem('fedgecolor',cgrey,'flonlimit',[-180 180]);
        gridm('gcolor',cgrey,'glinestyle','-');
        hold on
    end
    plot(xy_sm_left(:,1),xy_sm_left(:,2),'linewidth',3,'color',colorset(iF,:));
    
    subplot(122)  % right click
    if iF==1
        axesm eckert4
        axis off
        framem('fedgecolor',cgrey,'flonlimit',[-180 180]);
        gridm('gcolor',cgrey,'glinestyle','-');
        hold on
    end
    plot(xy_sm_right(:,1),xy_sm_right(:,2),'linewidth',3,'color',colorset(iF,:));
    
    clear xy_sm_*
end

subplot(121)
colormap(jet(num_freq))
colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
    'TickLabels',{num2str(param.freq_wanted'/1e3)},'location','southoutside');
grid
tightmap

subplot(122)
colormap(jet(num_freq))
colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
    'TickLabels',{num2str(param.freq_wanted'/1e3)},'location','southoutside');
grid
tightmap

% Save figures
save_fname = sprintf('%s_mf_cntr',data_path);
saveas(fig,fullfile(save_path,[save_fname,'.fig']),'fig');
saveSameSize(fig,'file',fullfile(save_path,[save_fname,'.png']),...
    'format','png','renderer','painters');

close(fig)

