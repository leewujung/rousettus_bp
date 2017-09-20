% 2015 12 28  Averaged elevation cut for individual clicks
% 2016 04 20  Update to work with new rotated data format
% 2016 05 07  Plot for NIFTI poster
% 2016 07 22  Plot for paper, need curves for different frequencies
%             Update loading rotate data part

clear

if isunix
    addpath('~/Dropbox/0_CODE/MATLAB/EllipseDirectFit');
    addpath('~/Dropbox/0_CODE/MATLAB/saveSameSize');
    
    data_base_path = '~/internal_2tb/Dropbox/0_ANALYSIS/bp_processing';
    results_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
else
    usrn = getenv('username');
    if strcmp(usrn,'Wu-Jung')
        addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
        addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    else
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
    end
    
    if strcmp(usrn,'Wu-Jung')
        data_base_path = 'F:\Dropbox\0_ANALYSIS\bp_processing\proc_output_rousettus_checked';
        results_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp\analysis_results_figs';
    else
        data_base_path = 'F:\Dropbox\0_ANALYSIS\bp_processing\proc_output_rousettus_checked';
        results_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp\analysis_results_figs';
    end
end

% el_cut_all = [1,3,15,30];
el_cut_all = 3;
freq_wanted_all = [20:5:60]*1e3;       % frequency to be plotted
save_opt = 1;
rotate_data_path = 'rotate_all_click_20160721';


% All rotate data files
rotate_file_all = dir(fullfile(results_base_path,rotate_data_path,'*.mat'));

for iE=1:length(el_cut_all)
    for iF=1:length(freq_wanted_all)
        freq_wanted = freq_wanted_all(iF);
        
        el_cut = el_cut_all(iE);  % elevation range to be averaged over
        if el_cut<3
            el_range = -el_cut:0.5:el_cut;
        else
            el_range = -el_cut:3:el_cut;
        end
        az_range = -180:3:180;
        [az_band,el_band] = meshgrid(az_range,el_range);
        tot_area = range(el_range)*range(az_range);
        
        % Set up various paths
        [~,script_name,~] = fileparts(mfilename('fullpath'));
        save_path = fullfile(results_base_path,script_name);
        if ~exist(save_path,'dir')
            mkdir(save_path);
        end
        
        cut_az = nan(length(rotate_file_all),size(el_band,2));
        click_side = nan(length(rotate_file_all),1);
        for iB = 1:length(rotate_file_all)
            % Load rotated data
            R = load(fullfile(results_base_path,...
                rotate_data_path,rotate_file_all(iB).name));
            
            % Load bp_proc data
            ss = strsplit(rotate_file_all(iB).name,'_');
            click_num = str2double(ss{7}(2:3));
            str = strjoin(ss(5:6),'_');
            bp_proc_file = sprintf('rousettus_20150825_%s_mic_data_bp_proc.mat',str);
            D = load(fullfile(data_base_path,bp_proc_file));
            
            % Get call_dB for specific frequency
            [call_dB,~,~,~] = get_call_azel_dB_data(D,freq_wanted,click_num);
            
            % Get other parameters
            call_dB_norm = call_dB-max(call_dB);
            az = R.shift_tilt_final.az;
            el = R.shift_tilt_final.el;
            ch_include_idx = (R.raw_meas.ch_include_idx(:) & ~isnan(az(:)))';
            vq_azel = rbfinterp([az_band(:)';el_band(:)'],...
                rbfcreate([az(ch_include_idx)';el(ch_include_idx)'],call_dB_norm(ch_include_idx),...
                'RBFFunction','multiquadrics'));
            vq_azel = reshape(vq_azel,size(az_band));
            
            % Set values outside of boundary to NaN
            azk = az(ch_include_idx);
            elk = el(ch_include_idx);
            k = boundary(azk,elk,0);  % outer boundary of all measured points
            [in,on] = inpolygon(az_band,el_band,azk(k),elk(k));
            in = in|on;
            vq_azel(~in) = NaN;
            
            % Get cut area
            az_band_cut = az_band(in);
            el_band_cut = el_band(in);
            kk = boundary(az_band_cut(:),el_band_cut(:),0);  % outer boundary of all measured points
            cut_area = polyarea(az_band_cut(kk),el_band_cut(kk));
            
            cut_az(iB,:) = nanmean(vq_azel,1);
            click_side(iB) = R.raw_meas.click_side;
            
        end  % loop through all clicks
        
        % Process left/right clicks separately
        cut_az_right = cut_az(click_side==1,:);
        cut_az_right_mean = nanmean(cut_az(click_side==1,:));
        cut_az_right_std = nanstd(cut_az(click_side==1,:));
        
        cut_az_left = cut_az(click_side==0,:);
        cut_az_left_mean = nanmean(cut_az(click_side==0,:));
        cut_az_left_std = nanstd(cut_az(click_side==0,:));
        
        % Plot
        % plot mean with all curves, separately
        fig_mean_all_curve = figure;
        subplot(211)
        plot(az_range,cut_az(click_side==1,:),'color',[1 1 1]*190/255)
        hold on
        hh = plot(az_range,cut_az_right_mean,'k','linewidth',1.5);
        ylim([-40 5])
        xlim([-180 180])
        grid
        ylabel('Relative intensity (dB)');
        ll = legend(hh,'Right');
        set(gca,'fontsize',12,'xtick',-180:60:180);
        subplot(212)
        plot(az_range,cut_az(click_side==0,:),'color',[1 1 1]*190/255)
        hold on
        hh = plot(az_range,cut_az_left_mean,'k','linewidth',1.5);
        ylim([-40 5])
        xlim([-180 180])
        grid
        ll = legend(hh,'Left');
        set(gca,'fontsize',12,'xtick',-180:60:180);
        xlabel('Azimuth (deg)');
        ylabel('Relative intensity (dB)');
        suptitle(sprintf('%s, el cut %d deg, freq %d kHz',...
            regexprep(script_name,'_','\\_'),el_cut,freq_wanted/1e3));
        if save_opt==1
            save_fname = sprintf('%s_elcut%02ddeg_%dkHz_mean_all_curves',...
                script_name,el_cut,freq_wanted/1e3);
            saveas(fig_mean_all_curve,fullfile(save_path,[save_fname,'.fig']),'fig');
            saveSameSize(fig_mean_all_curve,'file',fullfile(save_path,save_fname),...
                'format','png','renderer','painters');
        end
        close(fig_mean_all_curve)
        
    end  % loop through all freq
end  % loop through all el_cut
