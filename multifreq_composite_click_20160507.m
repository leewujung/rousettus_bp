% 2015 12 27  Make multi-freq composite clicks
% 2016 04 20  Update to make it compatible with new rotated click format
% 2016 05 07  Reprocess for plots for NIFTI poster

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
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_processing';
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
end

bat_proc_path = 'proc_output_rousettus_checked';

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_root = 'analysis_results_figs';
save_path = fullfile(base_path,save_root,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

data_path = 'rotate_all_click_20160419';

bat = {'3bat','36134','34271','39184'};

for iBAT=1%:length(bat)
data_file = ['rotate_all_click_20160419_all_click_rotated_',bat{iBAT},'.mat'];
D = load(fullfile(base_path,save_root,data_path,data_file));  % rotated data

A.param.base_path = base_path;
A.param.data_path = data_path;
A.param.data_file = data_file;
A.param.bp_proc_path = bat_proc_path;
A.param.rotated_data_path = data_path;
A.param.rotated_data_file = data_file;
A.param.map = D.map;

% Merge all data
num_ch = size(D.raw_meas.click_side_rep{1},2);
click_side_merge = reshape(cell2mat(D.raw_meas.click_side_rep),1,[]);
click_num_total = size(click_side_merge,2)/num_ch;  % total number of clicks included in the analysis
az_elpctr_tilt_all = reshape(cell2mat(D.shift_tilt_final.az),1,[]);  % az/el, rotated to ellipse center and compensated for rotation
el_elpctr_tilt_all = reshape(cell2mat(D.shift_tilt_final.el),1,[]);

A.param.num_ch = num_ch;
A.param.click_num_total = click_num_total;
A.scatter_data.all.click_side = click_side_merge;
A.scatter_data.all.az_shift_tilt = az_elpctr_tilt_all;
A.scatter_data.all.el_shift_tilt = el_elpctr_tilt_all;

% Get all freq data
% freq_wanted = (25:5:55)*1e3;
freq_wanted = (20:5:50)*1e3;
num_freq = length(freq_wanted);

path.base_path = base_path;
path.data_path = bat_proc_path;
multifreq_data = get_multifreq_data(path,D.bp_processed_file_all,freq_wanted);

A.param.freq_wanted = freq_wanted;
A.multifreq_data = multifreq_data;


% Composite click params
threshold = 3;  % only use averaged results if number of measurement > threshold
binsize = 10;  % bin size in azimuth and elevation

A.param.composite_threshold = threshold;
A.param.composite_binsize = binsize;


% right clicks
az_right = az_elpctr_tilt_all(click_side_merge==1);
el_right = el_elpctr_tilt_all(click_side_merge==1);
idx_nan_right = isnan(az_right);  % index of NaN data point
click_num_right = length(az_right)/num_ch;
az_right(idx_nan_right) = [];
el_right(idx_nan_right) = [];

% left clicks
az_left = az_elpctr_tilt_all(click_side_merge==0);
el_left = el_elpctr_tilt_all(click_side_merge==0);
idx_nan_left = isnan(az_left);  % index of NaN data point
click_num_left = length(az_left)/num_ch;
az_left(idx_nan_left) = [];
el_left(idx_nan_left) = [];

A.scatter_data.right.idx_nan = idx_nan_right;
A.scatter_data.right.az_shift_tilt = az_right;
A.scatter_data.right.el_shift_tilt = el_right;
A.scatter_data.left.idx_nan = idx_nan_left;
A.scatter_data.left.az_shift_tilt = az_left;
A.scatter_data.left.el_shift_tilt = el_left;

% initialize
call_dB_norm_merge = nan(num_freq,length(click_side_merge));
c3db_xy_all_right = cell(num_freq,1);
c3db_xy_all_left = cell(num_freq,1);
call_dB_right = zeros(num_freq,sum(~idx_nan_right));
call_dB_left = zeros(num_freq,sum(~idx_nan_left));
for iF=1:num_freq
    fprintf('Processing data at %d kHz\n',freq_wanted(iF)/1e3);
    
    % All data points for a single freq
    call_dB_norm_merge(iF,:) = reshape(cell2mat(multifreq_data.call_dB_norm(:,iF)),1,[]);
    
    tmp = call_dB_norm_merge(iF,click_side_merge==1);
    tmp(idx_nan_right) = [];  % take out NaN data point
    call_dB_right(iF,:) = tmp;
    tmp = call_dB_norm_merge(iF,click_side_merge==0);
    tmp(idx_nan_left) = [];  % take out NaN data point
    call_dB_left(iF,:) = tmp;
    
    % Interpolation using all data points
    [~,vq_norm_all_right(iF,:,:),azq_right,elq_right] = interp_bp(az_right/180*pi,el_right/180*pi,call_dB_right(iF,:),'rbf');
    [~,vq_norm_all_left(iF,:,:),azq_left,elq_left] = interp_bp(az_left/180*pi,el_left/180*pi,call_dB_left(iF,:),'rbf');
    
    % Interpolation using averaged data
    [interp_avg_right(iF),bin_avg_right(iF)] = average_call(az_right,el_right,call_dB_right(iF,:),binsize,'eckert4',threshold);
    [interp_avg_left(iF),bin_avg_left(iF)] = average_call(az_left,el_left,call_dB_left(iF,:),binsize,'eckert4',threshold);
    
    % Get -3dB contour from averaged data
    [azq,elq,vq] = get_ortho_grid_azel(interp_avg_left(iF).azq_avg,interp_avg_left(iF).elq_avg,interp_avg_left(iF).vq_norm_avg);
    [~,c_level_nan_left] = get_main_contour(vq,unique(azq(:)),unique(elq(:)),-3);
    [c3db_xy_left(:,1),c3db_xy_left(:,2)] = mfwdtran(D.map.mstruct,c_level_nan_left(:,2),c_level_nan_left(:,1));  % [az,el] to [x,y]
    
    [azq,elq,vq] = get_ortho_grid_azel(interp_avg_right(iF).azq_avg,interp_avg_right(iF).elq_avg,interp_avg_right(iF).vq_norm_avg);
    [~,c_level_nan_right] = get_main_contour(vq,unique(azq(:)),unique(elq(:)),-3);
    [c3db_xy_right(:,1),c3db_xy_right(:,2)] = mfwdtran(D.map.mstruct,c_level_nan_right(:,2),c_level_nan_right(:,1));  % [az,el] to [x,y]
    
    c3db_xy_all_right{iF} = c3db_xy_right;
    c3db_xy_all_left{iF} = c3db_xy_left;
    
    clear c3db_xy_left c3db_xy_right
end

A.scatter_data.all.call_dB_norm = call_dB_norm_merge;
A.scatter_data.right.call_dB_norm = call_dB_right;
A.scatter_data.left.call_dB_norm = call_dB_left;

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
    save_fname = sprintf('%s_all_clicks_%s.mat',script_name,bat{iBAT});
    save(fullfile(save_path,save_fname),'-struct','A');
end

clear vq_norm_all_* interp_avg_* bin_avg_* call_dB_* c3db_xy_all_*

end  % loop through all bats