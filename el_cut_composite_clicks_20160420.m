% 2016 04 20  Adapt el_cut_indiv_clicks.m for composite clicks

clear

usrn = getenv('username');
if strcmp(usrn,'Wu-Jung')
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
else
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
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

% el_cut_all = 15;
el_cut_all = [1,3,15,30];
% freq_wanted = 35e3;       % frequency to be plotted

for iE=1:length(el_cut_all)
    el_cut = el_cut_all(iE);  % elevation range to be averaged over
    if el_cut<3
        el_range = -el_cut:0.5:el_cut;
    else
        el_range = -el_cut:3:el_cut;
    end
    az_range = -180:3:180;
    [az_band,el_band] = meshgrid(az_range,el_range);
    tot_area = range(el_range)*range(az_range);
    
    % Each individual bat
    bat = {'3bat','36134','34271','39184'};
    for iBAT=2%1:length(bat)
        data_file = ['multifreq_composite_click_20160420_all_clicks_',bat{iBAT},'.mat'];
        D = load(fullfile(base_path,save_root,data_path,data_file));
        
        num_freq = length(D.param.freq_wanted);
        for iF=1:num_freq
            vq_azel_left = rbfinterp([az_band(:)';el_band(:)'],...
                rbfcreate([D.averaged_composite.left.bin(iF).az_avg(:)';...
                D.averaged_composite.left.bin(iF).el_avg(:)'],...
                D.averaged_composite.left.bin(iF).avg_call_dB(:)',...
                'RBFFunction','multiquadrics'));
            vq_azel_left = reshape(vq_azel_left,size(az_band));
            vq_azel_left = mean(vq_azel_left,1);
            
            vq_azel_right = rbfinterp([az_band(:)';el_band(:)'],...
                rbfcreate([D.averaged_composite.right.bin(iF).az_avg(:)';...
                D.averaged_composite.right.bin(iF).el_avg(:)'],...
                D.averaged_composite.right.bin(iF).avg_call_dB(:)',...
                'RBFFunction','multiquadrics'));
            vq_azel_right = reshape(vq_azel_right,size(az_band));
            vq_azel_right = mean(vq_azel_right,1);
            
            % Plot left and right composite curves separately
            fig_composite_el_cut = figure;
            hright = plot(az_range,vq_azel_right,'color','r','linewidth',1.5);
            hold on
            hleft = plot(az_range,vq_azel_left,'color','b','linewidth',1.5);
            ylim([-40 5])
            xlim([-180 180])
            grid
            ll = legend([hright,hleft],{'Right','Left'});
            set(ll,'fontsize',12);
            xlabel('Azimuth (deg)');
            ylabel('Relative intensity (dB)');
            title(sprintf('%s, %s, cut %d deg, freq %02d kHz',...
                regexprep(script_name,'_','\\_'),bat{iBAT},el_cut,D.param.freq_wanted(iF)/1e3));
            if save_plot_opt==1
                save_fname = sprintf('%s_elcut%02ddeg_%s_f%02dkHz.png',...
                    script_name,el_cut,bat{iBAT},D.param.freq_wanted(iF)/1e3);
                saveSameSize(fig_composite_el_cut,'file',fullfile(save_path,save_fname),...
                    'format','png','renderer','painters');
            end
            close(fig_composite_el_cut)
        end
    end  % loop through all bats
    
end  % loop through all el_cut
