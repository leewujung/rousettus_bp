% 2015 12 28  Plots for multi-freq composite clicks
% 2016 04 20  Separate this part out from multifreq_composite_click.m
% 2016 05 07  Plot for NIFTI poster

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
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_processing';
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
end

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_root = 'analysis_results_figs';
save_path = fullfile(base_path,save_root,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

data_path = 'multifreq_composite_click_20160420';

% Each individual bat
bat = {'3bat','36134','34271','39184'};
for iBAT=1%:length(bat)
    data_file = ['multifreq_composite_click_20160420_all_clicks_',bat{iBAT},'.mat'];
    load(fullfile(base_path,save_root,data_path,data_file));

    num_freq = length(param.freq_wanted);
    for iF=3%1:num_freq
        
        vq_norm_min = -27;
        contour_vec = 0:-3:(floor(vq_norm_min/3)-1)*3;
        cvec_min_idx = find(contour_vec-vq_norm_min<0,1,'first');
        
        fig_avg_bp = figure;
        set(fig_avg_bp,'position',[60,190,1250,390]);
        
        subplot(121)  % left click
        axesm eckert4
        framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
        gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
        axis off
        contourfm(averaged_composite.left.interp(iF).elq_avg,...
            averaged_composite.left.interp(iF).azq_avg,...
            averaged_composite.left.interp(iF).vq_norm_avg,...
            contour_vec(2:cvec_min_idx),...
            'fill','on','linecolor','w');  % don't plot 0 dB contour
%         contourfm(averaged_composite.left.interp(iF).elq_avg,...
%             averaged_composite.left.interp(iF).azq_avg,...
%             averaged_composite.left.interp(iF).vq_norm_avg,-30:3:-3);
        colorbar('southoutside','ticks',fliplr(contour_vec(1:cvec_min_idx)));
        colormap(parula(cvec_min_idx-1));
        caxis([contour_vec(cvec_min_idx) 0]);
        tightmap
%         caxis([-30 0])
%         colormap(parula(10))
%         colorbar('Ticks',-30:3:0,'location','southoutside');
        title(sprintf('Averaged left clicks, bat %s, %d kHz, th=%d, bin=%ddeg',...
            bat{iBAT},param.freq_wanted(iF)/1e3,param.composite_threshold,param.composite_binsize));
        
        subplot(122)  % right click
        axesm eckert4
        framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
        gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
        axis off
        contourfm(averaged_composite.right.interp(iF).elq_avg,...
            averaged_composite.right.interp(iF).azq_avg,...
            averaged_composite.right.interp(iF).vq_norm_avg,...
            contour_vec(2:cvec_min_idx),...
            'fill','on','linecolor','w');  % don't plot 0 dB contour
%         contourfm(averaged_composite.right.interp(iF).elq_avg,...
%             averaged_composite.right.interp(iF).azq_avg,...
%             averaged_composite.right.interp(iF).vq_norm_avg,-3:-3:-30);
        colorbar('southoutside','ticks',fliplr(contour_vec(1:cvec_min_idx)));
        colormap(parula(cvec_min_idx-1));
        caxis([contour_vec(cvec_min_idx) 0]);
        tightmap
%         caxis([-30 0])
%         colormap(parula(10))
%         colorbar('Ticks',-30:3:0,'location','southoutside');
        title(sprintf('Averaged right clicks, bat %s, %d kHz, th=%d, bin=%ddeg',...
            bat{iBAT},param.freq_wanted(iF)/1e3,param.composite_threshold,param.composite_binsize));
        
        if save_plot_opt == 1
            save_fname = sprintf('%s_%s_f%02dkHz',script_name,bat{iBAT},param.freq_wanted(iF)/1e3);
            saveas(fig_avg_bp,fullfile(save_path,[save_fname,'.fig']),'fig');
            saveSameSize_300(fig_avg_bp,'file',fullfile(save_path,save_fname),...
                'format','png','renderer','painters');
        end
        close(fig_avg_bp)
    end
end


