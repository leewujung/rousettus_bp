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
    for iK=1:knum
        % Mean curves
        subplot(knum+1,2,1);
        plot(C.az_range,nanmean(C.cut_az_left(param.idx_cluster_left==iK,:),1),'linewidth',2);
        hold on
        subplot(knum+1,2,2);
        plot(C.az_range,nanmean(C.cut_az_right(param.idx_cluster_right==iK,:),1),'linewidth',2);
        hold on
        
        % Left click
        subplot(knum+1,2,3+2*(iK-1))
        axesm eckert4
        framem('fedgecolor',cgrey,'flonlimit',[-180 180]);
        gridm('gcolor',cgrey,'glinestyle','-');
        axis off
%         surfm(averaged_composite.left.interp(iF,iK).elq_avg,...
%             averaged_composite.left.interp(iF,iK).azq_avg,...
%             averaged_composite.left.interp(iF,iK).vq_norm_avg);
        contourfm(averaged_composite.left.interp(iF,iK).elq_avg,...
            averaged_composite.left.interp(iF,iK).azq_avg,...
            averaged_composite.left.interp(iF,iK).vq_norm_avg,...
            contour_vec(2:cvec_min_idx),...
            'fill','on','linecolor','w');  % don't plot 0 dB contour
        colorbar('eastoutside','ticks',fliplr(contour_vec(1:cvec_min_idx)));
        colormap(parula(cvec_min_idx-1));
        caxis([contour_vec(cvec_min_idx) 0]);
        tightmap
        title(sprintf('Left, cluster %d',iK));
        
        % Right click
        subplot(knum+1,2,4+2*(iK-1));  
        axesm eckert4
        framem('fedgecolor',cgrey,'flonlimit',[-180 180]);
        gridm('gcolor',cgrey,'glinestyle','-');
        axis off
%         surfm(averaged_composite.right.interp(iF,iK).elq_avg,...
%             averaged_composite.right.interp(iF,iK).azq_avg,...
%             averaged_composite.right.interp(iF,iK).vq_norm_avg);
        contourfm(averaged_composite.right.interp(iF,iK).elq_avg,...
            averaged_composite.right.interp(iF,iK).azq_avg,...
            averaged_composite.right.interp(iF,iK).vq_norm_avg,...
            contour_vec(2:cvec_min_idx),...
            'fill','on','linecolor','w');  % don't plot 0 dB contour
        colorbar('eastoutside','ticks',fliplr(contour_vec(1:cvec_min_idx)));
        colormap(parula(cvec_min_idx-1));
        caxis([contour_vec(cvec_min_idx) 0]);
        tightmap
        title(sprintf('Right, cluster %d',iK));
        
    end % cluster loop
    
    subplot(knum+1,2,1);
    xlabel('Azimuth (deg)'); ylabel('Norm''ed energy');
    xlim([-180 180]); ylim([-30 3]); grid on
    legend({num2str([1:knum]')});
    subplot(knum+1,2,2);
    xlabel('Azimuth (deg)'); ylabel('Norm''ed energy');
    xlim([-180 180]); ylim([-30 3]); grid on
    legend({num2str([1:knum]')});
    
    suptitle(sprintf('Avg bp clustered at %dkHz, bat %s, %dkHz th=%d, bin=%d, knum=%d',...
        param.freq_used_for_cluster/1e3,bat_num,param.freq_wanted(iF)/1e3,...
        param.composite_threshold,param.composite_binsize,knum));

    % Save figures
    save_fname = sprintf('%s_knum%d_%dkHz',...
        data_file,knum,param.freq_wanted(iF)/1e3);
    %     saveas(fig,fullfile(save_path,[save_fname,'.fig']),'fig');
    saveSameSize(fig,'file',fullfile(save_path,[save_fname,'.png']),...
        'format','png','renderer','painters');
    
    close(fig)

% end % freq loop


