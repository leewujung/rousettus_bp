% 2015 11 24  Use real measurement az/el angles for model beampattern
% 2016 04 19  Use measured mean azimuth width to infer piston model diameter
%             Most processing and plotting codes are from rotate_all_click_20160419.m
% 2016 05 10  Re-process for plots in paper
% 2016 07 21  Add noise at different levels to model
% 2016 07 25  Set flag when shift_rotate_bp error out, see notes for detail
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

noise_std_all = [1,0,2];  % [dB]

for iN=1:length(noise_std_all)

rng_seed = 168;
noise_mean = 0;  % added noise profile [dB]
noise_std = noise_std_all(iN);
results_path = 'analysis_results_figs';
rotate_data_path = 'rotate_all_click_20161024';
match = 'azel';  % az|el|azel

rng(rng_seed);  % seed the random number generator

% Find all click files to match
rotate_data_file_all = dir(fullfile(data_base_path,results_path,rotate_data_path,'*.mat'));


% Set save path
[~,script_name,~] = fileparts(mfilename('fullpath'));
script_name = sprintf('%s_%s_std%2.1f',script_name,match,noise_std);
save_path = fullfile(data_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

azel_distr_path = 'az_el_distr_model_data_cmp_20160722';
azel_distr_file = 'az_el_distr_model_data_cmp_20160722_data.mat';
azel_distr = load(fullfile(data_base_path,results_path,azel_distr_path,azel_distr_file));

S.param.base_path = data_base_path;
S.param.rotate_data_path = rotate_data_path;
S.param.rotate_data_file_all = rotate_data_file_all;
S.param.azel_distr_path = azel_distr_path;
S.param.azel_distr_file = azel_distr_file;
S.param.noise_mean = noise_mean;
S.param.noise_std = noise_std;
S.param.rng_seed = rng_seed;
S.param.match = match;


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


% Set bp param
bp_info.c = 344;  % sound speed [m/s]
bp_info.freq = freq_model;  % [Hz]
bp_info.k = 2*pi*bp_info.freq/bp_info.c;  % wavenumber
bp_info.type = 'piston';


% Find aperture with best-fitting -3dB contour to mean azimuth of data
a_fine = (2:0.01:7)*1e-3; % [m]
theta_fine = 0:pi/1000:pi/2;
piston_bp = 20*log10(abs(2*besselj(1,bp_info.k*a_fine'*sin(theta_fine))./...
    (bp_info.k*a_fine'*sin(theta_fine))));
[C,~] = contour(theta_fine/pi*180,a_fine,piston_bp,[-3 -6],'fill','on');
c_level = parse_contour_output(C);
piston_3dB_c = [c_level(2).X; c_level(2).Y]';

switch match
    case 'az'
        match_deg = mean(azel_distr.az_data);
    case 'el'
        match_deg = mean(azel_distr.el_data);
    case 'azel'
        match_deg = mean([azel_distr.az_data;azel_distr.el_data]);
end

[~,idx_mean] = min(abs(piston_3dB_c(:,1)-match_deg));
a_mean = piston_3dB_c(idx_mean,2);
close
bp_info.a = a_mean;

S.model_beampattern_info = bp_info;


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
    az = data.raw.az/180*pi;  % convert to [rad]
    el = data.raw.el/180*pi;
    
    [~,mmidx] = max(call_dB);
    call_max_azel = [az(mmidx),el(mmidx)];
    
    % Model mic output
    piston_bp = model_beam_dB(bp_info,call_max_azel,[az el]);
    piston_bp = piston_bp+randn(size(piston_bp))*noise_std+noise_mean;
    
    % Fit ellipse
    [click_side,raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
        shift_rotate_bp(az,el,piston_bp','eckert4',it_shift_th);  % modified 2016/10/25
    
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
        title(regexprep(save_fname,'_','\\_'));
        saveSameSize(fig_fit,'file',...
            fullfile(save_path,sprintf('%s_contour.png',save_fname)),...
            'format','png','renderer','painters');
        close(fig_fit)
    end
    
    % Rotated model output
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
            suptitle(sprintf('%s, %s, all clicks',script_name, regexprep(t_name,'_','\\_')));
            suptitle([save_fname,' all clicks']);
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

