% 2015 12 27  Make multi-freq composite clicks

clear

bat='39184';
species='rousettus';

addpath('E:\repositories\beampattern_processing\fit_beams\');
addpath('..\')

save_opt = 1;

% Load compiled rotated data
base_path = '..\..\';

compile_path = 'proc_output_beam_align';
compile_file = ['all_click_rotated_' bat '.mat'];
load(fullfile(base_path,compile_path,compile_file));

% Path/filename for saving data and figures
save_header = bat;
save_path = fullfile(base_path,'composite_figs');  % set path for saving files
if ~exist(save_path,'dir')
  mkdir(save_path);
end

A.param.compiled_data_path = compile_path;
A.param.compiled_data_file = compile_file;
A.param.map.map_projection = map_projection;
A.param.map.map_plot_xlim = map_plot_xlim;
A.param.map.map_plot_ylim = map_plot_ylim;
A.param.map.mstruct = mstruct;


%% Merge all data
num_ch = size(raw_meas.click_side_rep{1},2);
click_side_merge = reshape(cell2mat(raw_meas.click_side_rep),1,[]);
click_num_total = size(click_side_merge,2)/num_ch;  % total number of clicks included in the analysis
az_elpctr_tilt_all = reshape(cell2mat(shift_tilt_final.az),1,[]);  % az/el, rotated to ellipse center and compensated for rotation
el_elpctr_tilt_all = reshape(cell2mat(shift_tilt_final.el),1,[]);

A.param.num_ch = num_ch;
A.param.click_num_total = click_num_total;
A.scatter_data.all.click_side = click_side_merge;
A.scatter_data.all.az_shift_tilt = az_elpctr_tilt_all;
A.scatter_data.all.el_shift_tilt = el_elpctr_tilt_all;

% Get all freq data
freq_wanted = (30:5:70)*1e3;
num_freq = length(freq_wanted);
path.base_path = base_path;
path.data_path = data_path;
multifreq_data = get_multifreq_data(path,processed_files,freq_wanted);

A.param.freq_wanted = freq_wanted;
A.multifreq_data = multifreq_data;


% Composite click params
threshold = 3;  % only use averaged results if number of measurement > threshold
binsize = 10;  % bin size in azimuth and elevation

A.param.composite_threshold = threshold;
A.param.composite_binsize = binsize;

if strcmp(species,'rousettus')
  % right clicks
  az_right = az_elpctr_tilt_all(click_side_merge==1);
  el_right = el_elpctr_tilt_all(click_side_merge==1);
  idx_nan_right = isnan(az_right);  % index of NaN data point
  click_num_right = size(az_right)/num_ch;
  az_right(idx_nan_right) = [];
  el_right(idx_nan_right) = [];
  
  % left clicks
  az_left = az_elpctr_tilt_all(click_side_merge==0);
  el_left = el_elpctr_tilt_all(click_side_merge==0);
  idx_nan_left = isnan(az_left);  % index of NaN data point
  click_num_left = size(az_left)/num_ch;
  az_left(idx_nan_left) = [];
  el_left(idx_nan_left) = [];
  
  A.scatter_data.right.idx_nan = idx_nan_right;
  A.scatter_data.right.az_shift_tilt = az_right;
  A.scatter_data.right.el_shift_tilt = el_right;
  A.scatter_data.left.idx_nan = idx_nan_left;
  A.scatter_data.left.az_shift_tilt = az_left;
  A.scatter_data.left.el_shift_tilt = el_left;
  
  % [radian] to [deg] conversion
  % az_right = az_right/pi*180;
  % el_right = el_right/pi*180;
  % az_left = az_left/pi*180;
  % el_left = el_left/pi*180;
  
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
    [c3db_xy_left(:,1),c3db_xy_left(:,2)] = mfwdtran(mstruct,c_level_nan_left(:,2),c_level_nan_left(:,1));  % [az,el] to [x,y]
    
    [azq,elq,vq] = get_ortho_grid_azel(interp_avg_right(iF).azq_avg,interp_avg_right(iF).elq_avg,interp_avg_right(iF).vq_norm_avg);
    [~,c_level_nan_right] = get_main_contour(vq,unique(azq(:)),unique(elq(:)),-3);
    [c3db_xy_right(:,1),c3db_xy_right(:,2)] = mfwdtran(mstruct,c_level_nan_right(:,2),c_level_nan_right(:,1));  % [az,el] to [x,y]
    
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
  
else
  %no left/right stuff
  
  az_all = az_elpctr_tilt_all(:);
  el_all = el_elpctr_tilt_all(:);
  idx_nan_all = isnan(az_all);  % index of NaN data point
  click_num_all = length(az_all)/num_ch;
  az_all(idx_nan_all) = [];
  el_all(idx_nan_all) = [];
  
  A.scatter_data.all.idx_nan = idx_nan_all;
  A.scatter_data.all.az_shift_tilt = az_all;
  A.scatter_data.all.el_shift_tilt = el_all;
  
  % [radian] to [deg] conversion
  % az_all = az_all/pi*180;
  % el_all = el_all/pi*180;
  
  % initialize
  call_dB_norm_merge = nan(num_freq,length(az_elpctr_tilt_all));
  c3db_xy_all_all = cell(num_freq,1);
  call_dB_all = zeros(num_freq,sum(~idx_nan_all));
  for iF=1:num_freq
    fprintf('Processing data at %d kHz\n',freq_wanted(iF)/1e3);
    
    % All data points for a single freq
    call_dB_norm_merge(iF,:) = reshape(cell2mat(multifreq_data.call_dB_norm(:,iF)),1,[]);
    
    tmp = call_dB_norm_merge(iF,:);
    tmp(idx_nan_all) = [];  % take out NaN data point
    call_dB_all(iF,:) = tmp;
        
    % Interpolation using all data points
    [~,vq_norm_all_all(iF,:,:),azq_all,elq_all] = ...
      interp_bp(az_all/180*pi,el_all/180*pi,call_dB_all(iF,:),'rbf');
    
    % Interpolation using averaged data
    [interp_avg_all(iF),bin_avg_all(iF)] = ...
      average_call(az_all,el_all,call_dB_all(iF,:),binsize,'eckert4',threshold);
    
    % Get -3dB contour from averaged data
    [azq,elq,vq] = get_ortho_grid_azel(interp_avg_all(iF).azq_avg,...
      interp_avg_all(iF).elq_avg,interp_avg_all(iF).vq_norm_avg);
    [~,c_level_nan_all] = get_main_contour(vq,unique(azq(:)),unique(elq(:)),-3);
    [c3db_xy_all(:,1),c3db_xy_all(:,2)] = ...
      mfwdtran(mstruct,c_level_nan_all(:,2),c_level_nan_all(:,1));  % [az,el] to [x,y]
    
    c3db_xy_all_all{iF} = c3db_xy_all;
    
    clear c3db_xy_all
  end
  
  A.scatter_data.all.call_dB_norm = call_dB_norm_merge;
  A.scatter_data.all.call_dB_norm = call_dB_all;
  
  A.raw_composite.all.vq = vq_norm_all_all;  % raw composite all click
  A.raw_composite.all.azq = azq_all;
  A.raw_composite.all.elq = elq_all;
  
  A.averaged_composite.all.interp = interp_avg_all;
  A.averaged_composite.all.bin = bin_avg_all;
  
  A.multifreq_3dB_contour.all = c3db_xy_all_all;
end


if save_opt==1
  save(fullfile(save_path,[save_header,'_composite_data.mat']),'-struct','A');
end


