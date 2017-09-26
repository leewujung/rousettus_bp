% 2017 09 25  Check the effects of missing head markers regarding estimation
%             of head aim 

clear
warning off
if isunix
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/rbfinterp_v1.2');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/EllipseDirectFit');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/saveSameSize');
    base_path = '~/internal_2tb/Dropbox/0_ANALYSIS/bp_processing';
    save_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
else
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\EllipseDirectFit');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_processing';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
end


% Set up various paths
[~,script_name,~] = fileparts(mfilename('fullpath'));
results_path = 'analysis_results_figs';
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

bat_proc_path = 'proc_output_rousettus_new_checked';
bat_proc_file_all = dir(fullfile(base_path,bat_proc_path,'rousettus_20150825_*'));

% Process all files
for iB = 1:length(bat_proc_file_all)

    bat_proc_file = bat_proc_file_all(iB).name;
    ss = strsplit(bat_proc_file,'_');
    save_fname = strjoin([script_name,ss(3:4)],'_');

    fprintf('File: %s\n',bat_proc_file);
    
    data = load(fullfile(base_path,bat_proc_path,bat_proc_file));
    S(iB).call_num = length(data.proc.chk_good_call);
    S(iB).good_call_idx = data.proc.chk_good_call';
    S(iB).source_head_aim = data.proc.source_head_aim;

    % load from processed bp data
    sm_len = data.track.smooth_len;
    track_t = data.track.track_raw_time;
    track_int_t = data.track.track_interp_time;
    track_sm = data.track.track_smooth;
    track_int = data.track.track_interp;
    marker_inc = data.track.marker_indicator;
    head_n_int = data.head_normal.head_normal_int;
    head_aim_int = data.head_aim.head_aim_int;

    % generate fake head normal from mics on the floor
    A = data.mic_loc(data.param.mic_floor_idx,:);
    A0 = bsxfun(@minus,A,mean(A,1)); % Subtract "mean" point
    [~,~,V] = svd(A0,0);
    % --- vector normalization
    norm_vec = norm_mtx_vec(V(:,3)');
    % --- fake head aim/normal from track
    head_aim_fake = [track_sm(sm_len:end,1:2)-track_sm(1:end-sm_len+1,1:2), ...
                     zeros(size(track_sm,1)-sm_len+1,1)];
    % --- vector normalization
    head_aim_fake = norm_mtx_vec(head_aim_fake);
    % --- interpolation
    %head_aim_fake_int = int_track(track_t(1:end-sm_len+1),head_aim_fake,track_int_t);
    head_aim_fake_int(:,1) = interp1(track_t(1:end-sm_len+1),head_aim_fake(:,1),track_int_t);
    head_aim_fake_int(:,2) = interp1(track_t(1:end-sm_len+1),head_aim_fake(:,2),track_int_t);
    head_aim_fake_int(:,3) = interp1(track_t(1:end-sm_len+1),head_aim_fake(:,3),track_int_t);

    % identify calls to be compared
    notnanidx = ~isnan(head_aim_int(:,1));  % these are the good calls to be compared
    head_aim_int_fake = nan(size(head_aim_int));
    head_n_int_fake = nan(size(head_n_int));
    head_aim_int_fake(notnanidx,:) = head_aim_fake_int(notnanidx,:);
    head_n_int_fake(notnanidx,:) = repmat(norm_vec,sum(notnanidx),1);

    % Get calls 
    proc_call_num = length(data.mic_data.call_idx_w_track);
    aim_v = nan(proc_call_num,3);
    aim_v_fake = nan(proc_call_num,3);
    norm_v = nan(proc_call_num,3);
    norm_v_fake = nan(proc_call_num,3);
    for iC=1:proc_call_num
        curr_call_loc_idx_on_track = data.track.call_loc_idx_on_track_interp(iC);
        aim_v(iC,:) = head_aim_int(curr_call_loc_idx_on_track,:);  % head aim at call
        norm_v(iC,:) = head_n_int(curr_call_loc_idx_on_track,:);  % head normal at call
        aim_v_fake(iC,:) = head_aim_int_fake(curr_call_loc_idx_on_track,:);  % head aim at call
        norm_v_fake(iC,:) = head_n_int_fake(curr_call_loc_idx_on_track,:);  % head normal at call
    end

    idx_comp = data.proc.chk_good_call==1 & data.proc.source_head_aim'==1;
    if any(idx_comp)
        aim_v = aim_v(idx_comp,:);
        norm_v = norm_v(idx_comp,:);
        aim_v_fake = aim_v_fake(idx_comp,:);
        norm_v_fake = norm_v_fake(idx_comp,:);

        % Calculate angular offset between real and fake aim and head normal vectors
        S(iB).aim_theta_deg_xyz = acosd(dot(aim_v',aim_v_fake'));  % each vector already normalized
        S(iB).norm_theta_deg_xyz = acosd(dot(norm_v',norm_v_fake'));  % each vector already normalized

        % X-Y plan vector direction only
        aim_v_2d = norm_mtx_vec(aim_v(:,1:2));
        aim_v_fake_2d = norm_mtx_vec(aim_v_fake(:,1:2));
        S(iB).aim_theta_deg_xy = acosd(dot(aim_v_2d',aim_v_fake_2d'));  % each vector already normalized
    end
end

% Plot quality check
fig_chk = figure('units','normalized','outerposition',[0 0 1 1]);
img_h = max([S(:).call_num]);
for iB = 1:length(bat_proc_file_all)
    img = zeros(img_h,2);
    img(1:S(iB).call_num,:) = 1;
    img(S(iB).good_call_idx==1,1) = 2;
    img(S(iB).source_head_aim==1,2) = 2;
    subplot(1,length(bat_proc_file_all),iB);
    imagesc(img);
    title(sprintf('%02d',iB),'fontsize',16);
    caxis([0,2])
end
hp = get(subplot(1,length(bat_proc_file_all),length(bat_proc_file_all)),'Position');
colormap(parula(3))
colorbar('Position', [0.92,0.2,0.01,0.65],...
         'ticks',[0.33,1,1.66],'ticklabels',{'No data','Bad','Good'},...
         'fontsize',16);
saveas(gcf,fullfile(save_path,[script_name,'_goodidx.png']),'png');
saveas(gcf,fullfile(save_path,[script_name,'_goodidx.fig']),'fig');



% Percentage of good calls with head marker info
all_good_call_idx = [S(:).good_call_idx];  % if good call: 1-yes, 0-no
all_source_head_aim = [S(:).source_head_aim]; % head aim from: 1-marker, 0-track

good_call_with_marker = sum(all_good_call_idx & all_source_head_aim);
good_call_no_marker = sum(all_good_call_idx & ~all_source_head_aim);

perc_with_marker = good_call_with_marker/(good_call_with_marker+good_call_no_marker);



% Stats of the deviation of head aim angles
aim_theta_deg_xyz_all = [S(:).aim_theta_deg_xyz];
norm_theta_deg_xyz_all = [S(:).norm_theta_deg_xyz];
aim_theta_deg_xy_all = [S(:).aim_theta_deg_xy];

fig_theta = figure;
corder = get(gca,'colororder');
subplot(211)
histogram(aim_theta_deg_xyz_all,0:5:100,'normalization','probability');
title('Head aim difference in 3D');
ylabel('Relative frequency');
subplot(212)
histogram(aim_theta_deg_xy_all,0:5:100,'normalization','probability');
ylim([0,0.2])
title('Head aim difference in 2D');
ylabel('Relative frequency');
xlabel('Degree');
saveas(gcf,fullfile(save_path,[script_name,'_aim_diff.png']),'png');
saveas(gcf,fullfile(save_path,[script_name,'_aim_diff.fig']),'fig');
