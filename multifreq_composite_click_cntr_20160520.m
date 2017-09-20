% 2015 12 28  Plots for multi-freq composite clicks
% 2016 04 20  Separate this part out from multifreq_composite_click.m
% 2016 05 07  Plot for NIFTI poster
% 2016 05 07  Plot for ASA talk

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

data_path = 'multifreq_composite_click_20160507';

% Each individual bat
fig_left = figure;
fig_right = figure;
bat = {'3bat','36134','34271','39184'};
for iBAT=1%:length(bat)
    data_file = ['multifreq_composite_click_20160507_all_clicks_',bat{iBAT},'.mat'];
    load(fullfile(base_path,save_root,data_path,data_file));
    
    num_freq = length(param.freq_wanted);
    colorset = jet(num_freq);
    contour_sm_len = 10;
    
    for iF=1:num_freq
        xy_sm_left(:,1) = smooth(multifreq_3dB_contour.left{iF}(:,1),contour_sm_len);
        xy_sm_left(:,2) = smooth(multifreq_3dB_contour.left{iF}(:,2),contour_sm_len);
        xy_sm_left(isnan(xy_sm_left(:,1)),:) = NaN;
        
        xy_sm_right(:,1) = smooth(multifreq_3dB_contour.right{iF}(:,1),contour_sm_len);
        xy_sm_right(:,2) = smooth(multifreq_3dB_contour.right{iF}(:,2),contour_sm_len);
        xy_sm_right(isnan(xy_sm_right(:,1)),:) = NaN;
        
        figure(fig_left)  % left click
        if iF==1
            axesm eckert4
            axis off
            framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
            gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
            hold on
        end
        plot(xy_sm_left(:,1),xy_sm_left(:,2),'linewidth',3,'color',colorset(iF,:));
        %         plot(multifreq_3dB_contour.left{iF}(:,1),multifreq_3dB_contour.left{iF}(:,2),'linewidth',2,'color',colorset(iF,:));
        figure(fig_right)  % right click
        if iF==1
            axesm eckert4
            axis off
            framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
            gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
            hold on
        end
        plot(xy_sm_right(:,1),xy_sm_right(:,2),'linewidth',3,'color',colorset(iF,:));
        %         plot(multifreq_3dB_contour.right{iF}(:,1),multifreq_3dB_contour.right{iF}(:,2),'linewidth',2,'color',colorset(iF,:));
        
        clear xy_sm_*
    end
    figure(fig_left)
    colormap(jet(num_freq))
    colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
        'TickLabels',{num2str(param.freq_wanted'/1e3)},'location','southoutside');
    grid
    tightmap
    title(sprintf('Multi-frequency contour, bat %s',bat{iBAT}));
    
    figure(fig_right)
    colormap(jet(num_freq))
    colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
        'TickLabels',{num2str(param.freq_wanted'/1e3)},'location','southoutside');
    grid
    tightmap
    title(sprintf('Multi-frequency contour, bat %s',bat{iBAT}));
    
    if save_plot_opt == 1
        save_fname_left = sprintf('%s_%s_left',script_name,bat{iBAT});
        saveas(fig_left,...
            fullfile(save_path,sprintf('%s.fig',save_fname_left)),'fig');
        saveSameSize(fig_left,'file',fullfile(save_path,save_fname_left),...
            'format','png','renderer','painters');
        
        save_fname_right = sprintf('%s_%s_rigth',script_name,bat{iBAT});
        saveas(fig_right,...
            fullfile(save_path,sprintf('%s.fig',save_fname_right)),'fig');
        saveSameSize(fig_right,'file',fullfile(save_path,save_fname_right),...
            'format','png','renderer','painters');
    end
end

