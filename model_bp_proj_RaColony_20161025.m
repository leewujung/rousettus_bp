% 2015 11 24  Use real measurement az/el angles for model beampattern
% 2016 04 19  Use measured mean azimuth width to infer piston model diameter
%             Most processing and plotting codes are from rotate_all_click_20160419.m
% 2016 05 10  Re-process for plots in paper
% 2016 07 21  Add noise at different levels to model
% 2016 07 25  Set flag when shift_rotate_bp error out, see notes for detail
% 2016 08 04  Adapt the code for projecting phased array model using
%             bullethead
% 2016 10 25  Update for version 1025 with out-of-bound points

clear

if isunix
    addpath('~/Dropbox/0_CODE/MATLAB/rbfinterp_v1.2');
    addpath('~/Dropbox/0_CODE/MATLAB/EllipseDirectFit');
    addpath('~/Dropbox/0_CODE/MATLAB/saveSameSize');
    addpath('~/Dropbox/0_CODE/beampattern_processing');
    
    data_base_path = '~/internal_2tb/Dropbox/0_ANALYSIS/bp_processing';
    save_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
else
    usrn = getenv('username');
    if strcmp(usrn,'Wu-Jung')
        addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
        addpath('F:\Dropbox\0_CODE\MATLAB\EllipseDirectFit');
        addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
        addpath('F:\Dropbox\0_CODE\beampattern_processing');
    else
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\EllipseDirectFit']);
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\beampattern_processing']);
    end
    
    % Set various paths
    if strcmp(usrn,'Wu-Jung')
        %         data_base_path = 'F:\Dropbox\0_ANALYSIS\bp_processing';
        data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
        save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    else
        data_base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
        save_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
    end
end

noise_std_all = [1,2];  % [dB]
results_path = 'analysis_results_figs';
rotate_data_path = 'rotate_all_click_20161024';
rng_seed = 168;

rng(rng_seed);  % seed the random number generator

% Param
freq_model = 35e3;
it_shift_th = 0.005;
azel_bnd = [300,150];  % bound of which the rotate-shift difference larger than this would be discarded
cvec = 0:-3:-39;

plot_opt_all_click = 1;
plot_opt_indiv_click = 1;
save_opt = 1;

S.param.freq_modeled = freq_model;
S.param.iterative_shift_threshold = it_shift_th;
S.param.shift_bound_azel = azel_bnd;


% Set map projection
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);

S.map.map_projection = map_proj;
S.map.mstruct = mstruct;

% Set bp model prediction path/files
bp_model_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp\bp_bem_modeling\model_bp_save_20161009';
bp_model_file = 'model_bp_save_20161009_Ra-colony-rotear-0.5mm_35kHz_x029_y000_z-05.mat';
BP = load(fullfile(bp_model_path,bp_model_file));

% Load model prediction
idxnotnan = ~isnan(BP.pp_plot);
[~,BP.vq_norm,BP.azq,BP.elq] = ...
    interp_bp(BP.az(idxnotnan)/180*pi,BP.el(idxnotnan)/180*pi,BP.pp_plot(idxnotnan),'rbf');
BP.azq = BP.azq/pi*180;
BP.elq = BP.elq/pi*180;
[~,mm_idx] = max(BP.vq_norm(:));
[BP.max_i,BP.max_j] = ind2sub(size(BP.vq_norm),mm_idx);

S.bp_model_path = bp_model_path;
S.bp_model_file = bp_model_file;


% Find all click files to match
rotate_data_file_all = dir(fullfile(data_base_path,results_path,rotate_data_path,'*.mat'));


for iN=1:length(noise_std_all)

rng_seed = 168;
noise_mean = 0;  % added noise profile [dB]
noise_std = noise_std_all(iN);


% Set save path
[~,script_name,~] = fileparts(mfilename('fullpath'));
script_name = sprintf('%s_std%2.1f',script_name,noise_std);
save_path = fullfile(data_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

S.param.base_path = data_base_path;
S.param.rotate_data_path = rotate_data_path;
S.param.rotate_data_file_all = rotate_data_file_all;
S.param.noise_mean = noise_mean;
S.param.noise_std = noise_std;
S.param.rng_seed = rng_seed;


% Model piston projection
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
    call_dB = data.raw.call_dB;
    az = data.raw.az;  % convert to [rad]
    el = data.raw.el;
    side = data.raw_meas.click_side;  % 1-right, 0-left
    
    [~,mmidx] = max(call_dB);
    az_max = az(mmidx);
    el_max = el(mmidx);

    [~,vq_norm,azq,elq] = interp_bp(az/180*pi,el/180*pi,call_dB,'rbf');
    azq = azq/pi*180;
    elq = elq/pi*180;
    
    % Model bp
    if side==1
        azq_model = BP.azq;
        az_model = BP.az;
    else
        azq_model = -BP.azq;
        az_model = -BP.az;
    end
    elq_model = BP.elq;
    el_model = BP.el;
    azq_model_max = azq_model(BP.max_i,BP.max_j);
    elq_model_max = elq_model(BP.max_i,BP.max_j);

    % Rotate model bp according to max az/el
    el_diff = el_max-elq_model_max;
    az_diff = az_max-azq_model_max;
    [elq_model_rot,azq_model_rot] = rotatem(elq_model,azq_model,...
        [el_diff,az_diff],'inverse','degrees');
    [el_model_rot,az_model_rot] = rotatem(el_model,az_model,...
        [el_diff,az_diff],'inverse','degrees');
    
    % Project model bp to mic loc
    idxnotnan = ~isnan(BP.pp_plot);
    vq_mic = rbfinterp([az(:)';el(:)'],...
        rbfcreate([az_model_rot(idxnotnan)';el_model_rot(idxnotnan)'],BP.pp_plot(idxnotnan)',...
        'RBFFunction','multiquadrics'));
    
    % Add noise
    vq_mic = vq_mic+randn(size(vq_mic))*noise_std+noise_mean;

    % Interpolate according to mic loc
    [~,vq_mic_norm,azq_mic,elq_mic] = interp_bp(az/180*pi,el/180*pi,vq_mic,'rbf');
    azq_mic = azq_mic/pi*180;
    elq_mic = elq_mic/pi*180;
    
    % Fit ellipse
    [click_side,raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
        shift_rotate_bp(az/180*pi,el/180*pi,vq_mic,'eckert4',it_shift_th);
    
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
    S.click_side = click_side;
    S.raw = raw;
    S.rot_max = rot_max;
    S.rot_elpctr = rot_elpctr;
    S.rot_elpctr_tilt = rot_elpctr_tilt;
    
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
            saveSameSize(fig_clicks,'file',...
                fullfile(save_path,sprintf('%s_%s_all_clicks.png',script_name,t_name)),...
                'format','png','renderer','painters');
            close(fig_clicks)
        end
    end

    if save_opt==1
        save(fullfile(save_path,[save_fname,'.mat']),'-struct','S');
    end

end

end  % loop through all noise levels

