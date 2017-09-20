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
bat = {'36134','34271','39184'};
all_az_shift_left = [];
all_el_shift_left = [];
all_call_dB_norm_left = [];
all_az_shift_right = [];
all_el_shift_right = [];
all_call_dB_norm_right = [];
for iBAT=1:length(bat)
    data_file = ['multifreq_composite_click_20160420_all_clicks_',bat{iBAT},'.mat'];
    load(fullfile(base_path,save_root,data_path,data_file));

    num_freq = length(param.freq_wanted);
    for iF=1:num_freq
        % Plot for each individual bat
        fig_scat = figure;
        set(fig_scat,'position',[70,100,1200,560]);
        
        subplot(121)  % left click
        scatter(scatter_data.left.az_shift_tilt,...
            scatter_data.left.el_shift_tilt,15,...
            scatter_data.left.call_dB_norm(iF,:),'fill')
        xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
        title(sprintf('All left click data points, bat %s, %d kHz',bat{iBAT},param.freq_wanted(iF)/1e3));
        colorbar('location','southoutside');
        grid
        
        subplot(122)  % right click
        scatter(scatter_data.right.az_shift_tilt,...
            scatter_data.right.el_shift_tilt,15,...
            scatter_data.right.call_dB_norm(iF,:),'fill')
        xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
        title(sprintf('All right click data points, bat %s, %d kHz',bat{iBAT},param.freq_wanted(iF)/1e3));
        colorbar('location','southoutside');
        grid
        
        if save_plot_opt == 1
            save_fname = sprintf('%s_%s_f%02dkHz',script_name,bat{iBAT},param.freq_wanted(iF)/1e3);
            saveSameSize(fig_scat,'file',fullfile(save_path,save_fname),...
                'format','png','renderer','painters');
        end
        close(fig_scat)
    end
    % Merge data for all bats
    all_az_shift_left = [all_az_shift_left,scatter_data.left.az_shift_tilt];
    all_el_shift_left = [all_el_shift_left,scatter_data.left.el_shift_tilt];
    all_call_dB_norm_left = [all_call_dB_norm_left,scatter_data.left.call_dB_norm];
    all_az_shift_right = [all_az_shift_right,scatter_data.left.az_shift_tilt];
    all_el_shift_right = [all_el_shift_right,scatter_data.left.el_shift_tilt];
    all_call_dB_norm_right = [all_call_dB_norm_right,scatter_data.left.call_dB_norm];
end

% All bats merged
for iF=1:num_freq
    
    fig_scat_all = figure;
    set(fig_scat_all,'position',[70,100,1200,560]);
    
    subplot(121)  % left click
    scatter(all_az_shift_left,all_el_shift_left,15,all_call_dB_norm_left(iF,:),'fill');
    xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
    title(sprintf('All left click data points, 3 bats, %d kHz',param.freq_wanted(iF)/1e3));
    colorbar('location','southoutside');
    grid
    
    subplot(122)  % right click
    scatter(all_az_shift_right,all_el_shift_right,15,all_call_dB_norm_right(iF,:),'fill');    
    xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
    title(sprintf('All right click data points, 3 bats, %d kHz',param.freq_wanted(iF)/1e3));
    colorbar('location','southoutside');
    grid
    
    if save_plot_opt == 1
        save_fname = sprintf('%s_3bat_f%02dkHz',script_name,param.freq_wanted(iF)/1e3);
        saveSameSize(fig_scat_all,'file',fullfile(save_path,save_fname),...
            'format','png','renderer','painters');
    end
    
    close(fig_scat_all)
end


