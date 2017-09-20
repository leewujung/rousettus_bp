% 2015 12 28  Averaged elevation cut for individual clicks
% 2016 04 20  Update to work with new rotated data format
% 2016 05 07  Plot for NIFTI poster
% 2016 07 22  Plot for paper
%             Update loading rotate data part

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
        data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
        save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    else
        data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
        save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    end
end

% el_cut_all = [1,3,15,30];
el_cut_all = 3;
save_opt = 1;
results_path = 'analysis_results_figs';
rotate_data_path = 'rotate_all_click_20160721';

% All rotate data files
bat_proc_file_all = dir(fullfile(data_base_path,results_path,rotate_data_path,'*.mat'));


for iE=1:length(el_cut_all)

freq_wanted = 35e3;       % frequency to be plotted
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
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end


for iB = 1:length(bat_proc_file_all)
    D = load(fullfile(data_base_path,results_path,...
            rotate_data_path,bat_proc_file_all(iB).name));
    call_dB = D.raw_meas.call_dB;
    call_dB_norm = call_dB-max(call_dB);
    az = D.shift_tilt_final.az;
    el = D.shift_tilt_final.el;
    ch_include_idx = D.raw_meas.ch_include_idx(:) & ~isnan(az(:))';
    vq_azel = rbfinterp([az_band(:)';el_band(:)'],...
        rbfcreate([az(ch_include_idx);el(ch_include_idx)],call_dB_norm(ch_include_idx),...
        'RBFFunction','multiquadrics'));
    vq_azel = reshape(vq_azel,size(az_band));
    
    
    % Get averaged elevation band
    cut_az = cell(length(D.rotate_data),1);
    click_side = cell(length(D.rotate_data),1);
    for iB=1:length(D.rotate_data)
        num_good_call = length(D.raw_meas.good_call_idx{iB});
        cut_az{iB} = nan(num_good_call,length(az_range));
        click_side{iB} = nan(num_good_call,1);
        for iC=1:num_good_call
            call_dB = D.raw_meas.call_dB{iB}(iC,:);
            call_dB_norm = call_dB-max(call_dB);
            az = D.shift_tilt_final.az{iB}(iC,:);
            el = D.shift_tilt_final.el{iB}(iC,:);
            ch_include_idx = D.raw_meas.ch_include_idx{iB}(iC,:) & ~isnan(az);
            vq_azel = rbfinterp([az_band(:)';el_band(:)'],...
                rbfcreate([az(ch_include_idx);el(ch_include_idx)],call_dB_norm(ch_include_idx),...
                'RBFFunction','multiquadrics'));
            vq_azel = reshape(vq_azel,size(az_band));
            
            % Set values outside of boundary to NaN
            azk = az(ch_include_idx);
            elk = el(ch_include_idx);
            k = boundary(azk',elk',0);  % outer boundary of all measured points
            [in,on] = inpolygon(az_band,el_band,azk(k),elk(k));
            in = in|on;
            vq_azel(~in) = NaN;
            
            % Get cut area
            az_band_cut = az_band(in);
            el_band_cut = el_band(in);
            kk = boundary(az_band_cut(:),el_band_cut(:),0);  % outer boundary of all measured points
            cut_area = polyarea(az_band_cut(kk),el_band_cut(kk));
            
            %         if cut_area/tot_area>0.5  % criteria to include data
            cut_az{iB}(iC,:) = nanmean(vq_azel,1);
            click_side{iB}(iC) = D.raw_meas.click_side{iB}(iC);
            %         end
            
        end
    end
    
    % Process left/right clicks separately
    click_side_all = cell2mat(click_side);
    cut_az_all = cell2mat(cut_az);
    
    cut_az_right = cut_az_all(click_side_all==1,:);
    cut_az_left = cut_az_all(click_side_all==0,:);
    
    cut_az_right_mean = nanmean(cut_az_all(click_side_all==1,:));
    cut_az_left_mean = nanmean(cut_az_all(click_side_all==0,:));
    
    cut_az_right_std = nanstd(cut_az_all(click_side_all==1,:));
    cut_az_left_std = nanstd(cut_az_all(click_side_all==0,:));
    
    % Plot
    % plot mean with all curves, separately
    fig_mean_all_curve = figure;
    subplot(211)
    plot(az_range,cut_az_all(click_side_all==1,:),'color',[1 1 1]*190/255)
    hold on
    hh = plot(az_range,cut_az_right_mean,'k','linewidth',1.5);
    ylim([-40 5])
    xlim([-180 180])
    grid
    ylabel('Relative intensity (dB)');
    ll = legend(hh,'Right');
    set(gca,'fontsize',12,'xtick',-180:60:180);
    subplot(212)
    plot(az_range,cut_az_all(click_side_all==0,:),'color',[1 1 1]*190/255)
    hold on
    hh = plot(az_range,cut_az_left_mean,'k','linewidth',1.5);
    ylim([-40 5])
    xlim([-180 180])
    grid
    ll = legend(hh,'Left');
    set(gca,'fontsize',12,'xtick',-180:60:180);
    xlabel('Azimuth (deg)');
    ylabel('Relative intensity (dB)');
    suptitle(sprintf('%s, bat %s, el cut %d deg',regexprep(script_name,'_','\\_'),bat{iBAT},el_cut));
    if save_opt==1
        save_fname = sprintf('%s_elcut%02ddeg_%s_mean_all_curves',script_name,el_cut,bat{iBAT});
        saveas(fig_mean_all_curve,fullfile(save_path,[save_fname,'.fig']),'fig');
        saveSameSize(fig_mean_all_curve,'file',fullfile(save_path,save_fname),...
            'format','png','renderer','painters');
    end
    close(fig_mean_all_curve)
    
    
    % Plot mean with std, separately
    notnanidx_right = ~isnan(cut_az_right_mean) & ~isnan(cut_az_right_std);
    notnanidx_left = ~isnan(cut_az_left_mean) & ~isnan(cut_az_left_std);
    
    fig_mean_std = figure;
    subplot(211)
    patch([az_range(notnanidx_right),fliplr(az_range(notnanidx_right))],...
        [cut_az_right_mean(notnanidx_right)-cut_az_right_std(notnanidx_right),...
        fliplr(cut_az_right_mean(notnanidx_right)+cut_az_right_std(notnanidx_right))],...
        'r','edgecolor','none','facealpha',0.2);
    hold on
    hright = plot(az_range,cut_az_right_mean,'r','linewidth',2);
    grid
    ylim([-40 5])
    xlim([-180 180])
    ylabel('Relative intensity (dB)');
    set(gca,'fontsize',12,'xtick',-180:60:180);
    ll = legend(hright,'Right');
    set(ll,'fontsize',12);
    subplot(212)
    patch([az_range(notnanidx_left),fliplr(az_range(notnanidx_left))],...
        [cut_az_left_mean(notnanidx_left)-cut_az_left_std(notnanidx_left),...
        fliplr(cut_az_left_mean(notnanidx_left)+cut_az_left_std(notnanidx_left))],...
        'b','edgecolor','none','facealpha',0.2);
    hold on
    hleft = plot(az_range,cut_az_left_mean,'b','linewidth',2);
    grid
    ylim([-40 5])
    xlim([-180 180])
    xlabel('Azimuth (deg)');
    ylabel('Relative intensity (dB)');
    set(gca,'fontsize',12,'xtick',-180:60:180);
    ll = legend(hleft,'Left');
    suptitle(sprintf('%s, bat %s, el cut %d deg',regexprep(script_name,'_','\\_'),bat{iBAT},el_cut));
    if save_opt==1
        save_fname = sprintf('%s_elcut%02ddeg_%s_mean_std',script_name,el_cut,bat{iBAT});
        saveas(fig_mean_std,fullfile(save_path,[save_fname,'.fig']),'fig');
        saveSameSize(fig_mean_std,'file',fullfile(save_path,save_fname),...
            'format','png','renderer','painters');
    end
    close(fig_mean_std)
    
    % Plot mean with std, overlapping
    fig_mean_std_overlap = figure;
    patch([az_range(notnanidx_right),fliplr(az_range(notnanidx_right))],...
        [cut_az_right_mean(notnanidx_right)-cut_az_right_std(notnanidx_right),...
        fliplr(cut_az_right_mean(notnanidx_right)+cut_az_right_std(notnanidx_right))],...
        'r','edgecolor','none','facealpha',0.2);
    hold on
    hright = plot(az_range,cut_az_right_mean,'r','linewidth',2);
    patch([az_range(notnanidx_left),fliplr(az_range(notnanidx_left))],...
        [cut_az_left_mean(notnanidx_left)-cut_az_left_std(notnanidx_left),...
        fliplr(cut_az_left_mean(notnanidx_left)+cut_az_left_std(notnanidx_left))],...
        'b','edgecolor','none','facealpha',0.2);
    hleft = plot(az_range,cut_az_left_mean,'b','linewidth',2);
    grid
    ylim([-40 5])
    xlim([-180 180])
    ll = legend([hright,hleft],'Right','Left');
    set(ll,'fontsize',12);
    xlabel('Azimuth (deg)');
    ylabel('Relative intensity (dB)');
    set(gca,'fontsize',12,'xtick',-180:60:180);
    suptitle(sprintf('%s, bat %s, el cut %d deg',regexprep(script_name,'_','\\_'),bat{iBAT},el_cut));
    if save_opt==1
        save_fname = sprintf('%s_elcut%02ddeg_%s_mean_std_overlap',script_name,el_cut,bat{iBAT});
        saveas(fig_mean_std_overlap,fullfile(save_path,[save_fname,'.fig']),'fig');
        saveSameSize(fig_mean_std_overlap,'file',fullfile(save_path,save_fname),...
            'format','png','renderer','painters');
    end
    close(fig_mean_std_overlap)
    
end  % loop through all bats

end  % loop through all el_cut
