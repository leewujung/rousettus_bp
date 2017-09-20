% 2015 12 28  Averaged elevation cut for individual clicks
% 2016 04 20  Update to work with new rotated data format
% 2016 05 07  Plot for NIFTI poster
% 2016 07 22  Plot for paper, need curves for different frequencies
%             Update loading rotate data part
% 2016 08 09  Plot for paper; revised to save processed results
% 2016 08 18  Revise for plot for paper
% 2016 10 25  Update for version 1025 with out-of-bound points

clear

if isunix
    addpath('~/Dropbox/0_CODE/MATLAB/EllipseDirectFit');
    addpath('~/Dropbox/0_CODE/MATLAB/saveSameSize');
    
    data_base_path = '~/internal_2tb/Dropbox/0_ANALYSIS/bp_processing';
    save_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
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
        save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    else
        data_base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\proc_output_rousettus_checked'];
        save_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
    end
end

el_cut_all = 3;
freq_wanted_all = [25:10:55]*1e3;       % frequency to be plotted
rotate_data_path = 'rotate_all_click_20161024';
bat_num = 'all';  % all,34271,36134,39184

[~,script_name,~] = fileparts(mfilename('fullpath'));
results_path = 'analysis_results_figs';
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end


% All rotate data files
if strcmp(bat_num,'all')
    rotate_file_all = dir(fullfile(save_base_path,results_path,rotate_data_path,'*.mat'));
else
    rotate_file_all = dir(fullfile(save_base_path,results_path,rotate_data_path,['*',bat_num,'*.mat']));
end

A.rotate_file_all = rotate_file_all;

for iE=1:length(el_cut_all)
    
    for iF=1:length(freq_wanted_all)
        freq_wanted = freq_wanted_all(iF);
        
        A.freq = freq_wanted;
        
        el_cut = el_cut_all(iE);  % elevation range to be averaged over
        if el_cut<3
            el_range = -el_cut:0.5:el_cut;
        else
            el_range = -el_cut:3:el_cut;
        end
        az_range = -180:3:180;
        [az_band,el_band] = meshgrid(az_range,el_range);
        %tot_area = range(el_range)*range(az_range);

        A.el_cut = el_cut;
        A.el_range = el_range;
        A.az_range = az_range;
        A.az_band = az_band;
        A.el_band = el_band;

        cut_az = nan(length(rotate_file_all),size(el_band,2));
        click_side = nan(length(rotate_file_all),1);
        for iB = 1:length(rotate_file_all)
            % Load rotated data
            R = load(fullfile(save_base_path,results_path,...
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
        
        A.cut_az = cut_az;
        A.click_side = click_side;
        
        % Process left/right clicks separately
        cut_az_right = cut_az(click_side==1,:);
        cut_az_right_mean = nanmean(cut_az(click_side==1,:));
        cut_az_right_std = nanstd(cut_az(click_side==1,:));
        
        cut_az_left = cut_az(click_side==0,:);
        cut_az_left_mean = nanmean(cut_az(click_side==0,:));
        cut_az_left_std = nanstd(cut_az(click_side==0,:));
        
        A.cut_az_right = cut_az_right;
        A.cut_az_right_mean = cut_az_right_mean;
        A.cut_az_right_std = cut_az_right_std;
        A.cut_az_left = cut_az_left;
        A.cut_az_left_mean = cut_az_left_mean;
        A.cut_az_left_std = cut_az_left_std;
        
        % Plot
        % plot mean with all curves, separately
        fig_mean_all_curve = figure('position',[280 300 1400 350]);

        subplot(121)  % Left clicks
        plot(az_range,cut_az(click_side==0,:),'color',[1 1 1]*190/255)
        hold on
        hh = plot(az_range,cut_az_left_mean,'k','linewidth',1.5);
        ylim([-30 5])
        xlim([-180 180])
        grid
        ll = legend(hh,'Left');
        set(gca,'fontsize',12,'xtick',-180:60:180);
        xlabel('Azimuth (deg)');
        ylabel('Relative intensity (dB)');
 
        subplot(122)  % Right clicks
        plot(az_range,cut_az(click_side==1,:),'color',[1 1 1]*190/255)
        hold on
        hh = plot(az_range,cut_az_right_mean,'k','linewidth',1.5);
        ylim([-30 5])
        xlim([-180 180])
        grid
        xlabel('Azimuth (deg)');
        ll = legend(hh,'Right');
        set(gca,'fontsize',12,'xtick',-180:60:180);

        pause(0.05)
        
        mtit(sprintf('%s, el cut %d deg, freq %d kHz',...
            regexprep(script_name,'_','\\_'),el_cut,freq_wanted/1e3));
        
        save_fname = sprintf('%s_bat%s_elcut%02ddeg_%dkHz',...
            script_name,bat_num,el_cut,freq_wanted/1e3);
        save(fullfile(save_path,[save_fname,'.mat']),'-struct','A');
        saveas(fig_mean_all_curve,fullfile(save_path,[save_fname,'.fig']),'fig');
        saveSameSize(fig_mean_all_curve,'file',fullfile(save_path,save_fname),...
            'format','png','renderer','painters');
        epswrite(fullfile(save_path,[save_fname,'.eps']));
        close(fig_mean_all_curve)
        
    end  % loop through all freq
end  % loop through all el_cut
