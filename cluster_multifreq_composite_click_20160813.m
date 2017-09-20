% 2015 12 27  Make multi-freq composite clicks
% 2016 04 20  Update to make it compatible with new rotated click format
% 2016 05 07  Reprocess for plots for NIFTI poster
% 2016 08 08  Update rotated data format after 20160721
% 2016 08 13  Composite click according to clustering results

clear

usrn = getenv('username');
if strcmp(usrn,'Wu-Jung')
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
else
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
end

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    bat_proc_path = 'F:\Dropbox\0_ANALYSIS\bp_processing\proc_output_rousettus_checked';
    data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
else
    bat_proc_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\proc_output_rousettus_checked'];
    data_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
    save_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
end

results_path = 'analysis_results_figs';

% Composite click params
freq_wanted = (20:5:60)*1e3;
threshold = 0;  % only use averaged results if number of measurement > threshold
binsize = 10;  % bin size in azimuth and elevation
num_freq = length(freq_wanted);
bat_num = '36134';  % all,34271,36134,39184

A.param.freq_wanted = freq_wanted;
A.param.composite_threshold = threshold;
A.param.composite_binsize = binsize;

% Get data files
data_path = 'rotate_all_click_20160721';
if strcmp(bat_num,'all')
    data_file = dir(fullfile(data_base_path,results_path,data_path,'*.mat'));
else
    data_file = dir(fullfile(data_base_path,results_path,data_path,['*',bat_num,'*.mat']));
end

A.param.data_base_path = data_base_path;
A.param.rotated_data_path = data_path;
A.param.rotated_data_file = data_file;

% Load 1 file to get params
D = load(fullfile(data_base_path,results_path,data_path,data_file(1).name));  % rotated data
num_ch = length(D.raw_meas.call_dB);

A.param.num_ch = num_ch;
A.param.map = D.map;

% Set save folder
[~,script_name,~] = fileparts(mfilename('fullpath'));
script_name = sprintf('%s_bat%s_bin%d_th%d',script_name,bat_num,binsize,threshold);
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Clustering according to curve shape
knum = 3;
el_cut = 3;
freq_cluster = 35e3;
elcut_path = 'fig_el_cut_indiv_clicks_20160809';
elcut_file = sprintf('%s_bat%s_elcut%02ddeg_%02dkHz.mat',...
    elcut_path,bat_num,el_cut,freq_cluster/1e3);
C = load(fullfile(data_base_path,results_path,elcut_path,elcut_file));

cluster_az_idx_range = 41:81;
idx_cluster_left = kmeans(C.cut_az_left(:,cluster_az_idx_range),knum);
idx_cluster_right = kmeans(C.cut_az_right(:,cluster_az_idx_range),knum);

save_fname = sprintf('%s_knum%d_elcut%02d_freqcls%02dkHz',...
    script_name,knum,el_cut,freq_cluster/1e3);

A.param.elcut_path = elcut_path;
A.param.elcut_file = elcut_file;
A.param.knum = knum;
A.param.el_cut = el_cut;
A.param.freq_used_for_cluster = freq_cluster;
A.param.cluster_az_idx_range = cluster_az_idx_range;
A.param.idx_cluster_left = idx_cluster_left;
A.param.idx_cluster_right = idx_cluster_right;


% Compile data from all clicks
click_side_all = nan(length(data_file),num_ch);
az_elpctr_tilt_all = nan(length(data_file),num_ch);
el_elpctr_tilt_all = nan(length(data_file),num_ch);
call_dB_all = nan(num_freq,length(data_file),num_ch);
az_max_all = nan(length(data_file),num_ch);
el_max_all = nan(length(data_file),num_ch);
azq_max_all = nan(length(data_file),num_ch);
elq_max_all = nan(length(data_file),num_ch);
for iS=1:length(data_file)
    disp(['Loading ',data_file(iS).name]);
    D = load(fullfile(data_base_path,results_path,data_path,data_file(iS).name));  % rotated data
    
    % Load multi-freq data
    ss = strsplit(strtok(data_file(iS).name,'.'),'_');
    bat_num = ss{5};
    trial_num = ss{6};
    click_num = str2double(ss{7}(2:3));
    bat_proc_fname = sprintf('rousettus_20150825_%s_%s_mic_data_bp_proc.mat',bat_num,trial_num);
    data = load(fullfile(bat_proc_path,bat_proc_fname));
    for iF=1:num_freq
        [call_dB,~,~,ch_include_idx] = get_call_azel_dB_data(data,freq_wanted(iF),click_num);  % get good mic index
        call_dB = call_dB - max(call_dB);
        call_dB(~ch_include_idx) = NaN;
        call_dB_all(iF,iS,:) = call_dB;
    end
    
    % Max mic az/el location
    [~,mmidx] = max(D.raw.call_dB);
    [~,mmqidx] = max(D.raw.vq_norm(:));
    
    % Merge all data
    click_side_all(iS,:) = D.raw_meas.click_side*ones(1,length(num_ch));
    az_elpctr_tilt_all(iS,:) = D.rot_elpctr_tilt.az;
    el_elpctr_tilt_all(iS,:) = D.rot_elpctr_tilt.el;
    az_max_all(iS,1) = D.raw.az(mmidx);
    el_max_all(iS,1) = D.raw.el(mmidx);
    azq_max_all(iS,1) = D.raw.azq(mmqidx);
    elq_max_all(iS,1) = D.raw.elq(mmqidx);
end

A.scatter_data.click_side_all = click_side_all;
A.scatter_data.az_elpctr_tilt_all = az_elpctr_tilt_all;
A.scatter_data.el_elpctr_tilt_all = el_elpctr_tilt_all;
A.scatter_data.call_dB_all = call_dB_all;
A.scatter_data.az_max_all = az_max_all;
A.scatter_data.el_max_all = el_max_all;
A.scatter_data.azq_max_all = azq_max_all;
A.scatter_data.elq_max_all = elq_max_all;


% Composite bp across all frequencies
c3db_xy_all_right = cell(num_freq,knum);
c3db_xy_all_left = cell(num_freq,knum);
for iF=1:num_freq
    disp(['Processing data at ',num2str(freq_wanted(iF)/1e3),' kHz']);

    call_dB_merge = squeeze(call_dB_all(iF,:,:));
    idx_right = find(click_side_all==1);
    idx_left = find(click_side_all==0);
    for iK=1:knum
        idx_right_k = idx_right(idx_cluster_right==iK);
        idx_left_k = idx_left(idx_cluster_left==iK);
        
        call_dB_right_merge = reshape(call_dB_merge(idx_right_k,:),[],1);
        az_elpctr_tilt_right_merge = reshape(az_elpctr_tilt_all(idx_right_k,:),[],1);
        el_elpctr_tilt_right_merge = reshape(el_elpctr_tilt_all(idx_right_k,:),[],1);
        idx_good_right = ~isnan(az_elpctr_tilt_right_merge);
        call_dB_right = call_dB_right_merge(idx_good_right);
        az_right = az_elpctr_tilt_right_merge(idx_good_right);
        el_right = el_elpctr_tilt_right_merge(idx_good_right);
        
        call_dB_left_merge = reshape(call_dB_merge(idx_left_k,:),[],1);
        az_elpctr_tilt_left_merge = reshape(az_elpctr_tilt_all(idx_left_k,:),[],1);
        el_elpctr_tilt_left_merge = reshape(el_elpctr_tilt_all(idx_left_k,:),[],1);
        idx_good_left = ~isnan(az_elpctr_tilt_left_merge);
        call_dB_left = call_dB_left_merge(idx_good_left);
        az_left = az_elpctr_tilt_left_merge(idx_good_left);
        el_left = el_elpctr_tilt_left_merge(idx_good_left);
        
%         % Interpolation using all data points
%         [~,vq_norm_all_right(iF,iK,:,:),azq_right,elq_right] = ...
%             interp_bp(az_right'/180*pi,el_right'/180*pi,call_dB_right','rbf');
%         [~,vq_norm_all_left(iF,iK,:,:),azq_left,elq_left] = ...
%             interp_bp(az_left'/180*pi,el_left'/180*pi,call_dB_left','rbf');
        
        % Interpolation using averaged data
        [interp_avg_right(iF,iK),bin_avg_right(iF,iK)] = ...
            average_call(az_right,el_right,call_dB_right,binsize,'eckert4',threshold);
        [interp_avg_left(iF,iK),bin_avg_left(iF,iK)] = ...
            average_call(az_left,el_left,call_dB_left,binsize,'eckert4',threshold);
        
        % Get -3dB contour from averaged data
        [azq,elq,vq] = get_ortho_grid_azel(interp_avg_right(iF,iK).azq_avg,...
            interp_avg_right(iF,iK).elq_avg,interp_avg_right(iF,iK).vq_norm_avg);
        [~,c_level_nan_right] = get_main_contour(vq,unique(azq(:)),unique(elq(:)),-3);
        [c3db_xy_right(:,1),c3db_xy_right(:,2)] = ...
            mfwdtran(D.map.mstruct,c_level_nan_right(:,2),c_level_nan_right(:,1));  % [az,el] to [x,y]
        
        [azq,elq,vq] = get_ortho_grid_azel(interp_avg_left(iF,iK).azq_avg,...
            interp_avg_left(iF,iK).elq_avg,interp_avg_left(iF,iK).vq_norm_avg);
        [~,c_level_nan_left] = get_main_contour(vq,unique(azq(:)),unique(elq(:)),-3);
        [c3db_xy_left(:,1),c3db_xy_left(:,2)] = ...
            mfwdtran(D.map.mstruct,c_level_nan_left(:,2),c_level_nan_left(:,1));  % [az,el] to [x,y]

        c3db_xy_all_right{iF,iK} = c3db_xy_right;
        c3db_xy_all_left{iF,iK} = c3db_xy_left;
        
        clear c3db_xy_left c3db_xy_right

    end % cluster loop
        
end
                                
% A.raw_composite.right.vq = vq_norm_all_right;  % raw composite right click
% A.raw_composite.right.azq = azq_right;
% A.raw_composite.right.elq = elq_right;
% A.raw_composite.left.vq = vq_norm_all_left;  % raw composite left click
% A.raw_composite.left.azq = azq_left;
% A.raw_composite.left.elq = elq_left;

A.averaged_composite.right.interp = interp_avg_right;
A.averaged_composite.right.bin = bin_avg_right;
A.averaged_composite.left.interp = interp_avg_left;
A.averaged_composite.left.bin = bin_avg_left;

A.multifreq_3dB_contour.right = c3db_xy_all_right;
A.multifreq_3dB_contour.left = c3db_xy_all_left;

save(fullfile(save_path,save_fname),'-struct','A');


