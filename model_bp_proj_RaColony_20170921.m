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

noise_std_all = [1];  % [dB]
results_path = 'analysis_results_figs';
rotate_data_path = 'rotate_all_click_2010308';  % rotation on 20170308, note
                                                % the missing '7' in filename
rng_seed = 168;
bpctr_opt = 'ectr';  % 'max' -- max beam energy location
                    % 'top' -- averaged loc for all normalized beam energy>-1
                    % 'ectr' -- center of best-fitting ellipse

rng(rng_seed);  % seed the random number generator

% Set params
it_shift_th = 0.005;
cvec = 0:-3:-39;

S.param.iterative_shift_threshold = it_shift_th;

% Set map projection
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);

S.map.map_projection = map_proj;
S.map.mstruct = mstruct;

% Other params
plot_opt_all_click = 1;
plot_opt_indiv_click = 1;
save_opt = 1;

freq_model = 35e3;  % project 35 kHz model and also use 35 kHz for rotation
S.param.freq_modeled = freq_model;

% Find all click files to match
rotate_data_file_all = dir(fullfile(data_base_path,results_path,rotate_data_path,'*.mat'));

% Save params
S.param.rotate_data_path = rotate_data_path;
S.param.rotate_data_file_all = rotate_data_file_all;
S.param.rng_seed = rng_seed;


% Set bp model prediction path/files
bp_model_path = 'model_bp_save_20161009_multifreq';
bp_model_file = sprintf(['model_bp_save_20161009_multifreq_Ra-colony-rotear-' ...
                    '0.5mm_%dkHz_x029_y000_z-05.mat'],freq_model/1e3);
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
max_elq_loc = BP.elq(mm_idx);
max_azq_loc = BP.azq(mm_idx);
% --- averaged location of normalized beam energy >-1
xx(isnan(xx)) = -inf;
[~,sort_idx] = sort(xx,'descend');
ii = xx(sort_idx)>-1;
top_elq_loc = mean(BP.elq(sort_idx(ii)));
top_azq_loc = mean(BP.azq(sort_idx(ii)));
% --- fit ellipse
% ------- BP.el = bem_results.theta
% ------- BP.az = bem_results.phi
[raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
    shift_rotate_bp_composite(BP.az,BP.el,BP.pp,map_proj,0.005);
[el_ectr,az_ectr] = minvtran(mstruct,rot_max.E.x0,rot_max.E.y0);  % inverse map projection
[el_ectr_r,az_ectr_r] = rotatem(el_ectr,az_ectr,...
                                [max_elq_loc,max_azq_loc],...
                                'inverse','degrees');
BP.max_elq_loc = max_elq_loc;
BP.max_azq_loc = max_azq_loc;
BP.top_elq_loc = top_elq_loc;
BP.top_azq_loc = top_azq_loc;
BP.el_ectr_r = el_ectr_r;
BP.az_ectr_r = az_ectr_r;
S.BP = BP;


%% Loop through all noise level
for iN=1:length(noise_std_all)

    noise_mean = 0;  % added noise profile [dB]
    noise_std = noise_std_all(iN);

    % Set save path
    [~,script_name,~] = fileparts(mfilename('fullpath'));
    script_name = sprintf('%s_std%2.1f',script_name,noise_std);
    save_path = fullfile(data_base_path,results_path,script_name);
    if ~exist(save_path,'dir')
        mkdir(save_path);
    end

    S.param.noise_mean = noise_mean;
    S.param.noise_std = noise_std;

    % Model projection setup
    rotate_data = cell(length(rotate_data_file_all),1);
    flag_trial = 0;
    flag_count = 0;
    for iB = 1:length(rotate_data_file_all)
        
        bat_proc_file = rotate_data_file_all(iB).name;
        ss = strsplit(bat_proc_file,'_');
        save_fname = strjoin([script_name,ss(5:7)],'_');
        
        if ~flag_trial
            t_name = strjoin(ss(5:6),'_');
            trial_file_all = dir(fullfile(data_base_path,...
                                          results_path,rotate_data_path,['*_',t_name,'_*.mat']));
            flag_trial = 1;
        end
        
        if flag_count<length(trial_file_all)
            flag_count = flag_count+1;
        end
        
        if plot_opt_all_click && flag_count==1
            numrow = ceil(length(trial_file_all)/4);
            fig_clicks = figure;
            set(fig_clicks,'Position',[100 100 1400 240*numrow]);
        end
        
        fprintf('File: %s\n',bat_proc_file);
        
        % Get az/el from measurement
        data = load(fullfile(data_base_path,results_path,rotate_data_path,bat_proc_file));
        
        % Beam center from measurements
        % --- max beam energy point
        xx = data.raw.vq_norm(:);
        [~,mm_idx] = max(xx);
        max_elq_loc = data.raw.elq(mm_idx);  % [deg]
        max_azq_loc = data.raw.azq(mm_idx);
        % --- averaged location of point >-1 dB normalized beam energy
        xx(isnan(xx)) = -inf;
        [~,sort_idx] = sort(xx,'descend');
        ii = xx(sort_idx)>-1;
        top_elq_loc = mean(data.raw.elq(sort_idx(ii)));  % [deg]
        top_azq_loc = mean(data.raw.azq(sort_idx(ii)));
        % --- fit ellipse
        if isempty(data.rot_elpctr)  % if data.rot_max.E.x0-y0 are not within
                                     % valid minvtran range
            S.data = [];
            S.model_rot = [];

        else
            [el_ectr,az_ectr] = minvtran(mstruct,data.rot_max.E.x0,...
                                         data.rot_max.E.y0);  % inverse map projection
            [el_ectr_r,az_ectr_r] = rotatem(el_ectr,az_ectr,...  % [deg]
                                            [max_elq_loc,max_azq_loc],...
                                            'inverse','degrees');
            
            % Model bp
            if data.raw_meas.click_side==1  % right click --> no need to flip az
                azq_model = BP.azq;
                az_model = BP.az;
                max_azq_loc_model = BP.max_azq_loc;
                top_azq_loc_model = BP.top_azq_loc;
                az_ectr_r_model = BP.az_ectr_r;
            else   % left click
                azq_model = -BP.azq;
                az_model = -BP.az;
                max_azq_loc_model = -BP.max_azq_loc;
                top_azq_loc_model = -BP.top_azq_loc;
                az_ectr_r_model = -BP.az_ectr_r;
            end
            elq_model = BP.elq;
            el_model = BP.el;
            max_elq_loc_model = max_elq_loc;
            top_elq_loc_model = top_elq_loc;
            el_ectr_r_model = el_ectr_r;

            switch bpctr_opt
              case 'max'
                el_diff = max_elq_loc-max_elq_loc_model;
                az_diff = max_azq_loc-max_azq_loc_model;
              case 'top'
                el_diff = top_elq_loc-top_elq_loc_model;
                az_diff = top_azq_loc-top_azq_loc_model;
              case 'ectr'
                el_diff = el_ectr_r-el_ectr_r_model;
                az_diff = az_ectr_r-az_ectr_r_model;
            end
            
            S.data.max_elq_loc = max_elq_loc;
            S.data.max_azq_loc = max_azq_loc;
            S.data.top_elq_loc = top_elq_loc;
            S.data.top_azq_loc = top_azq_loc;
            S.data.el_ectr_r = el_ectr_r;
            S.data.az_ectr_r = az_ectr_r;
            S.data.el_diff = el_diff;
            S.data.az_diff = az_diff;
            
            % Rotate model bp according to max az/el
            [elq_model_rot,azq_model_rot] = rotatem(elq_model,azq_model,...
                                                    [el_diff,az_diff],'inverse','degrees');
            [el_model_rot,az_model_rot] = rotatem(el_model,az_model,...
                                                  [el_diff,az_diff],'inverse','degrees');
            
            % Project model bp to mic loc
            idxnotnan = ~isnan(BP.pp_plot);
            vq_mic = rbfinterp([data.raw.az(:)';data.raw.el(:)'],...
                               rbfcreate([az_model_rot(idxnotnan)';el_model_rot(idxnotnan)'],...
                                         BP.pp_plot(idxnotnan)',...
                                         'RBFFunction','multiquadrics'));

            % Add noise
            vq_mic = vq_mic+randn(size(vq_mic))*noise_std+noise_mean;

            % Fit ellipse
            [click_side,raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
                shift_rotate_bp(data.raw.az/180*pi,data.raw.el/180*pi,...
                                vq_mic,'eckert4',it_shift_th);
            
            % if the best-fitting ellipse for shift_max is outside of globe
            % then don't plot
            if isempty(rot_elpctr) || isempty(rot_elpctr_tilt)
                flag_elps_fail = 1;
            else
                flag_elps_fail = 0;
            end
            
            % Check az, el shift and tilt
            if plot_opt_indiv_click==1 && flag_elps_fail~=1
                Vpass.raw = raw;
                Vpass.rot_max = rot_max;
                Vpass.rot_elpctr = rot_elpctr;
                Vpass.rot_elpctr_tilt = rot_elpctr_tilt;
                
                % Movement of [az,el] of each mic
                fig_azel = plot_indiv_click_azel_movement(Vpass);
                title(regexprep(save_fname,'_','\\_'));
                saveSameSize(fig_azel,'file',...
                             fullfile(save_path,sprintf('%s_azel_chk.png',save_fname)),...
                             'format','png','renderer','painters');
                close(fig_azel)
                
                % Plot the shift/tilt procedure
                fig_fit = plot_indiv_click_rotate(Vpass,cvec,mstruct);
                suptitle(regexprep(save_fname,'_','\\_'));
                saveSameSize(fig_fit,'file',...
                             fullfile(save_path,sprintf('%s_contour.png',save_fname)),...
                             'format','png','renderer','painters');
                close(fig_fit)
            end
            
            % Rotated model output
            S.model_rot.click_side = click_side;
            S.model_rot.raw = raw;
            S.model_rot.rot_max = rot_max;
            S.model_rot.rot_elpctr = rot_elpctr;
            S.model_rot.rot_elpctr_tilt = rot_elpctr_tilt;
            
            % Plot rotated bp using ellipse center
            if flag_count<=length(trial_file_all) && flag_elps_fail~=1
                if plot_opt_all_click
                    figure(fig_clicks)
                    subplot(numrow,4,flag_count);
                    axesm(mstruct)
                    contour(rot_elpctr_tilt.xq,rot_elpctr_tilt.yq,rot_elpctr_tilt.vq_norm,cvec,'fill','on');
                    framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
                    gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
                    tightmap
                    axis off
                    title(sprintf('Call #%02d',str2double(ss{7}(2:3))));
                    colormap(parula(length(cvec)-1));
                    caxis([cvec(end), 0])
                end
            end

            if flag_count==length(trial_file_all)
                flag_trial = 0;
                flag_count = 0;
                if plot_opt_all_click
                    figure(fig_clicks)
                    suptitle(regexprep(sprintf('%s, %s, all clicks',script_name,t_name),'_','\\_'));
                    saveSameSize(fig_clicks,'file',...
                                 fullfile(save_path,sprintf('all_clicks_%s_%s.png',script_name,t_name)),...
                                 'format','png','renderer','painters');
                    close(fig_clicks)
                end
            end

        end  % test if data.rot_max.E.x0-y0 are in the valid range for minvtran

        if save_opt==1
            save(fullfile(save_path,[save_fname,'.mat']),'-struct','S');
        end

    end  % loop through all clicks

end  % loop through all noise levels


