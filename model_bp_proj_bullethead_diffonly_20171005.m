% 2015 11 24  Use real measurement az/el angles for model beampattern
% 2016 04 19  Use measured mean azimuth width to infer piston model diameter
%             Most processing and plotting codes are from rotate_all_click_20160419.m
% 2016 05 10  Re-process for plots in paper
% 2016 07 21  Add noise at different levels to model
% 2016 07 25  Set flag when shift_rotate_bp error out, see notes for detail
% 2016 08 04  Adapt the code for projecting phased array model using
%             bullethead
% 2016 10 25  Update for version 1025 with out-of-bound points
% 2017 09 21  -- Enable using different beam center point for projection
%             -- Saving info that will allow projection using multiple
%                frequencies in model_composite_RaColony3456_20170921.m
% 2017 09 22  Take out the part to actually carry out simulation
%             and only obtain critical params for making it happen.
%             The actual simulation will happen in model_composite_RaColony3456_20170921.m
%             This modification is done to enable multi-freq simulation
% 2017 10 05  Use the code from `model_bp_proj_RaColony_diffonly_20170921.m`
%             to do bullethead multi-freq simulation`

clear

if isunix
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/rbfinterp_v1.2');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/EllipseDirectFit');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/saveSameSize');
    addpath('~/internal_2tb/Dropbox/0_CODE/beampattern_processing');
    
    data_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
    save_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
    model_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp/bp_bem_modeling/';
else
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\EllipseDirectFit');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    addpath('F:\Dropbox\0_CODE\beampattern_processing');
    
    data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    model_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp\bp_bem_modeling';
end

% Set map projection
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);

S.map.map_projection = map_proj;
S.map.mstruct = mstruct;


% Set params
results_path = 'analysis_results_figs';
rotate_data_path = 'rotate_all_click_2010308';  % rotation on 20170308, note
                                                % the missing '7' in filename
rng_seed = 168;
bpctr_opt = 'ectr';  % 'max' -- max beam energy location
                     % 'top' -- averaged loc for all normalized beam energy>-1
                     % 'ectr' -- center of best-fitting ellipse
it_shift_th = 0.005;
cvec = 0:-3:-39;
freq_model = 35e3;  % project 35 kHz model and also use 35 kHz for rotation

% Find all click files to match
rotate_data_file_all = dir(fullfile(data_base_path,results_path,rotate_data_path,'*.mat'));

% Save params
S.param.data_base_path = data_base_path;
S.param.save_base_path = save_base_path;
S.param.model_base_path = model_base_path;
S.param.bpctr_opt = bpctr_opt;
S.param.freq_model = freq_model;
S.param.rotate_data_path = rotate_data_path;
S.param.rng_seed = rng_seed;
S.param.iterative_shift_threshold = it_shift_th;


% Other params
plot_opt_all_click = 1;
plot_opt_indiv_click = 1;
save_opt = 1;


% Set save path
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(data_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end


% Set bp model prediction path/files
bp_model_path = 'model_bp_save_20161013_multifreq';
bp_model_file = sprintf(['model_bp_save_20161013_multifreq_bullethead-sc-' ...
                    'colony-0.5mm_%dkHz_x028_y000_z-03.mat'],freq_model/1e3);
BP = load(fullfile(model_base_path,bp_model_path,bp_model_file));

BP.bp_model_path = bp_model_path;
BP.bp_model_file = bp_model_file;


% Load model prediction
idxnotnan = ~isnan(BP.pp_plot);
[~,BP.vq_norm,BP.azq,BP.elq] = ...
    interp_bp(BP.az(idxnotnan)/180*pi,BP.el(idxnotnan)/180*pi,BP.pp_plot(idxnotnan),'rbf');
BP.azq = BP.azq/pi*180;
BP.elq = BP.elq/pi*180;

% Find beam center for model
% --- max beam energy location
xx = BP.vq_norm(:);
[~,mm_idx] = max(xx);
BP.max.el = BP.elq(mm_idx);
BP.max.az = BP.azq(mm_idx);
% --- averaged location of normalized beam energy >-1
xx(isnan(xx)) = -inf;
[~,sort_idx] = sort(xx,'descend');
ii = xx(sort_idx)>-1;
BP.top.el = mean(BP.elq(sort_idx(ii)));
BP.top.az = mean(BP.azq(sort_idx(ii)));
% --- fit ellipse
% ------- BP.el = bem_results.theta
% ------- BP.az = bem_results.phi
[raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
    shift_rotate_bp_composite(BP.az,BP.el,BP.pp,map_proj,0.005);
[el_ectr,az_ectr] = minvtran(mstruct,rot_max.E.x0,rot_max.E.y0);  % inverse map projection
[BP.ectr.el,BP.ectr.az] = rotatem(el_ectr,az_ectr,...
                                  [BP.max.el,BP.max.az],...
                                  'inverse','degrees');
S.BP = BP;


% Model projection setup
rotate_data = cell(length(rotate_data_file_all),1);
for iB = 1:length(rotate_data_file_all)
    
    bat_proc_file = rotate_data_file_all(iB).name;
    ss = strsplit(bat_proc_file,'_');
    save_fname = strjoin([script_name,ss(5:7)],'_');
    
    fprintf('file %03d: %s\n',iB,bat_proc_file);

    S.param.rotate_data_file = bat_proc_file;
    
    % Get az/el from measurement
    data = load(fullfile(data_base_path,results_path,rotate_data_path, ...
                         bat_proc_file));
    S.raw_meas_from_mic = data.raw_meas;
    
    % Beam center from measurements
    % --- max beam energy point
    xx = data.raw.vq_norm(:);
    [~,mm_idx] = max(xx);
    meas.max.el = data.raw.elq(mm_idx);  % [deg]
    meas.max.az = data.raw.azq(mm_idx);
    % --- averaged location of point >-1 dB normalized beam energy
    xx(isnan(xx)) = -inf;
    [~,sort_idx] = sort(xx,'descend');
    ii = xx(sort_idx)>-1;
    meas.top.el = mean(data.raw.elq(sort_idx(ii)));  % [deg]
    meas.top.az = mean(data.raw.azq(sort_idx(ii)));
    % --- fit ellipse
    if isempty(data.rot_elpctr)  % if data.rot_max.E.x0-y0 are not within
                                 % valid minvtran range
        S.meas = [];
        S.diff = [];

    else
        [el_ectr,az_ectr] = minvtran(mstruct,data.rot_max.E.x0,...
                                     data.rot_max.E.y0);  % inverse map projection
        [meas.ectr.el,meas.ectr.az] = rotatem(el_ectr,az_ectr,...  % [deg]
                                              [meas.max.el,meas.max.az],...
                                              'inverse','degrees');

        % Determine beam centr
        switch bpctr_opt
          case 'max'
            bpctr_data = meas.max;
            bpctr_model = BP.max;
          case 'top'
            bpctr_data = meas.top;
            bpctr_model = BP.top;
          case 'ectr'
            bpctr_data = meas.ectr;
            bpctr_model = BP.ectr;
        end

        % Model bp
        if data.raw_meas.click_side==1  % right click --> no need to flip az
            az_model = BP.az;
            diff.az = bpctr_data.az - bpctr_model.az;
            diff.el = bpctr_data.el - bpctr_model.el;
        else   % left click --> need to flip az, and therefore becomes
               % '+' to get diff.az and diff.el
            az_model = -BP.az;
            diff.az = bpctr_data.az + bpctr_model.az;
            diff.el = bpctr_data.el + bpctr_model.el;
        end
        el_model = BP.el;

        S.meas = meas;
        S.diff = diff;
    end
    
    if save_opt==1
        save(fullfile(save_path,[save_fname,'.mat']),'-struct','S');
    end
    
end  % loop through all clicks



