% 2015 12 27  Make multi-freq composite clicks
% 2016 04 20  Update to make it compatible with new rotated click format
% 2016 05 07  Reprocess for plots for NIFTI poster
% 2016 08 08  Update rotated data format after 20160721
% 2016 10 25  Update for rotated clicks 20161025 which did not take out
%             out-of-bound points
% 2016 10 27  Revert to shift_rotate_bp_20161025a which doesn't
%             re-interpolate at every step

clear

usrn = getenv('username');
if strcmp(usrn,'Wu-Jung')
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
else
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
end

save_opt = 1;

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
threshold = 3;  % only use averaged results if number of measurement > threshold
binsize = 10;  % bin size in azimuth and elevation
num_freq = length(freq_wanted);
bat_num = 'all';  % all,34271,36134,39184

A.param.freq_wanted = freq_wanted;
A.param.composite_threshold = threshold;
A.param.composite_binsize = binsize;

% Get data files
data_path = 'rotate_all_click_20161027';
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

% Compile data from all clicks
click_side_all = nan(length(data_file),num_ch);
az_elpctr_tilt_all = nan(length(data_file),num_ch);
el_elpctr_tilt_all = nan(length(data_file),num_ch);
call_dB_all = nan(num_freq,length(data_file),num_ch);
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
    
    % Merge all data
    click_side_all(iS,:) = D.raw_meas.click_side*ones(1,length(num_ch));
    az_elpctr_tilt_all(iS,:) = D.rot_elpctr_tilt.az;
    el_elpctr_tilt_all(iS,:) = D.rot_elpctr_tilt.el;
end

A.scatter_data.click_side_all = click_side_all;
A.scatter_data.az_elpctr_tilt_all = az_elpctr_tilt_all;
A.scatter_data.el_elpctr_tilt_all = el_elpctr_tilt_all;
A.scatter_data.call_dB_all = call_dB_all;


% Composite bp across all frequencies
c3db_xy_all_right = cell(num_freq,1);
c3db_xy_all_left = cell(num_freq,1);
for iF=1:num_freq
    disp(['Processing data at ',num2str(freq_wanted(iF)/1e3),' kHz']);
    
    % All data points for a single freq
    call_dB_merge = reshape(squeeze(call_dB_all(iF,:,:)),[],1);
    click_side_merge = reshape(click_side_all,[],1);
    az_elpctr_tilt_merge = reshape(az_elpctr_tilt_all,[],1);
    el_elpctr_tilt_merge = reshape(el_elpctr_tilt_all,[],1);
    
    idx_right_good = find(click_side_merge==1 & ~isnan(call_dB_merge) &...
                          ~isnan(az_elpctr_tilt_merge));
    call_dB_right = call_dB_merge(idx_right_good);
    az_right = az_elpctr_tilt_all(idx_right_good);
    el_right = el_elpctr_tilt_all(idx_right_good);
    
    idx_left_good = find(click_side_merge==0 & ~isnan(call_dB_merge) &...
                         ~isnan(az_elpctr_tilt_merge));
    call_dB_left = call_dB_merge(idx_left_good);
    az_left = az_elpctr_tilt_all(idx_left_good);
    el_left = el_elpctr_tilt_all(idx_left_good);
    
    % Interpolation using all data points
    [~,vq_norm_all_right(iF,:,:),azq_right,elq_right] = ...
        interp_bp(az_right'/180*pi,el_right'/180*pi,call_dB_right','rbf');
    [~,vq_norm_all_left(iF,:,:),azq_left,elq_left] = ...
        interp_bp(az_left'/180*pi,el_left'/180*pi,call_dB_left','rbf');
    
    % Interpolation using averaged data
    [interp_avg_right(iF),bin_avg_right(iF)] = ...
        average_call(az_right,el_right,call_dB_right,binsize,'eckert4',threshold);
    [interp_avg_left(iF),bin_avg_left(iF)] = ...
        average_call(az_left,el_left,call_dB_left,binsize,'eckert4',threshold);
    
    % Get -3dB contour from averaged data
    [azq,elq,vq] = get_ortho_grid_azel(interp_avg_left(iF).azq_avg,...
        interp_avg_left(iF).elq_avg,interp_avg_left(iF).vq_norm_avg);
    [~,c_level_nan_left] = get_main_contour(vq,unique(azq(:)),unique(elq(:)),-3);
    [c3db_xy_left(:,1),c3db_xy_left(:,2)] = ...
        mfwdtran(D.map.mstruct,c_level_nan_left(:,2),c_level_nan_left(:,1));  % [az,el] to [x,y]
    
    [azq,elq,vq] = get_ortho_grid_azel(interp_avg_right(iF).azq_avg,...
        interp_avg_right(iF).elq_avg,interp_avg_right(iF).vq_norm_avg);
    [~,c_level_nan_right] = get_main_contour(vq,unique(azq(:)),unique(elq(:)),-3);
    [c3db_xy_right(:,1),c3db_xy_right(:,2)] = ...
        mfwdtran(D.map.mstruct,c_level_nan_right(:,2),c_level_nan_right(:,1));  % [az,el] to [x,y]
    
    c3db_xy_all_right{iF} = c3db_xy_right;
    c3db_xy_all_left{iF} = c3db_xy_left;
    
    clear c3db_xy_left c3db_xy_right
end
                                
A.raw_composite.right.vq = vq_norm_all_right;  % raw composite right click
A.raw_composite.right.azq = azq_right;
A.raw_composite.right.elq = elq_right;
A.raw_composite.left.vq = vq_norm_all_left;  % raw composite left click
A.raw_composite.left.azq = azq_left;
A.raw_composite.left.elq = elq_left;

A.averaged_composite.right.interp = interp_avg_right;
A.averaged_composite.right.bin = bin_avg_right;
A.averaged_composite.left.interp = interp_avg_left;
A.averaged_composite.left.bin = bin_avg_left;

A.multifreq_3dB_contour.right = c3db_xy_all_right;
A.multifreq_3dB_contour.left = c3db_xy_all_left;

if save_opt==1
    save_fname = [script_name,'_merged_clicks.mat'];
    save(fullfile(save_path,save_fname),'-struct','A');
end

