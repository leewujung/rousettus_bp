% 2015 12 28  Plots for multi-freq composite clicks
% 2016 04 20  Separate this part out from multifreq_composite_click.m

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
for iBAT=1:length(bat)
    data_file = ['multifreq_composite_click_20160420_all_clicks_',bat{iBAT},'.mat'];
    load(fullfile(base_path,save_root,data_path,data_file));
    
    num_freq = length(param.freq_wanted);
    colorset = jet(num_freq);
    contour_sm_len = 10;
    
    fig_multifreq_contour = figure;
    set(fig_multifreq_contour,'position',[150 150 1200 400]);
    for iF=1:num_freq
        xy_sm_left(:,1) = smooth(multifreq_3dB_contour.left{iF}(:,1),contour_sm_len);
        xy_sm_left(:,2) = smooth(multifreq_3dB_contour.left{iF}(:,2),contour_sm_len);
        xy_sm_left(isnan(xy_sm_left(:,1)),:) = NaN;
        
        xy_sm_right(:,1) = smooth(multifreq_3dB_contour.right{iF}(:,1),contour_sm_len);
        xy_sm_right(:,2) = smooth(multifreq_3dB_contour.right{iF}(:,2),contour_sm_len);
        xy_sm_right(isnan(xy_sm_right(:,1)),:) = NaN;
        
        subplot(121)  % left click
        if iF==1
            axesm eckert4
            axis off
            framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
            gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
            hold on
        end
        plot(xy_sm_left(:,1),xy_sm_left(:,2),'linewidth',2,'color',colorset(iF,:));
        %         plot(multifreq_3dB_contour.left{iF}(:,1),multifreq_3dB_contour.left{iF}(:,2),'linewidth',2,'color',colorset(iF,:));
        subplot(122)  % right click
        if iF==1
            axesm eckert4
            axis off
            framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
            gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
            hold on
        end
        plot(xy_sm_right(:,1),xy_sm_right(:,2),'linewidth',2,'color',colorset(iF,:));
        %         plot(multifreq_3dB_contour.right{iF}(:,1),multifreq_3dB_contour.right{iF}(:,2),'linewidth',2,'color',colorset(iF,:));
        
        clear xy_sm_*
    end
    subplot(121)
    colormap(jet(num_freq))
    colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
        'TickLabels',{num2str(param.freq_wanted'/1e3)},'location','southoutside');
    grid
    subplot(122)
    colormap(jet(num_freq))
    colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
        'TickLabels',{num2str(param.freq_wanted'/1e3)},'location','southoutside');
    grid
    suptitle(sprintf('Multi-frequency contour, bat %s',bat{iBAT}));
    
    if save_plot_opt == 1
        save_fname = sprintf('%s_%s',script_name,bat{iBAT});
        saveSameSize(fig_multifreq_contour,'file',fullfile(save_path,save_fname),...
            'format','png','renderer','painters');
    end
    close(fig_multifreq_contour)
end

