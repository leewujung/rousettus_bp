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
for iBAT=1%:length(bat)
    data_file = ['multifreq_composite_click_20160420_all_clicks_',bat{iBAT},'.mat'];
    load(fullfile(base_path,save_root,data_path,data_file));

    num_freq = length(param.freq_wanted);
    for iF=1:num_freq
        
        fig_comp_raw = figure;
        set(fig_comp_raw,'position',[60,190,1250,390]);
        
        subplot(121)  % left click
        axesm eckert4
        framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
        gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
        axis off
        contourfm(raw_composite.left.elq/pi*180,...
            raw_composite.left.azq/pi*180,...
            squeeze(raw_composite.left.vq(iF,:,:)),-30:3:-3);
        caxis([-30 0])
        colormap(parula(10))
        colorbar('Ticks',-30:3:0,'location','southoutside');
        title(sprintf('Merged left clicks, bat %s, %d kHz, no averaging',bat{iBAT},param.freq_wanted(iF)/1e3));
        
        subplot(122)  % right click
        axesm eckert4
        framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
        gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
        axis off
        contourfm(raw_composite.right.elq/pi*180,...
            raw_composite.right.azq/pi*180,...
            squeeze(raw_composite.right.vq(iF,:,:)),-3:-3:-30);
        caxis([-30 0])
        colormap(parula(10))
        colorbar('Ticks',-30:3:0,'location','southoutside');
        title(sprintf('Merged right clicks, bat %s, %d kHz, no averaging',bat{iBAT},param.freq_wanted(iF)/1e3));
        
        if save_plot_opt == 1
            save_fname = sprintf('%s_%s_f%02dkHz',script_name,bat{iBAT},param.freq_wanted(iF)/1e3);
            saveSameSize(fig_comp_raw,'file',fullfile(save_path,save_fname),...
                'format','png','renderer','painters');
        end
        close(fig_comp_raw)
    end
end


