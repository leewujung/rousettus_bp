% 2016 08 08  Assemble composite clicks from simulated mic receptions
% 2016 10 25  Update for version 1025 with out-of-bound points
% 2017 09 22  Modififed from `model_composite_RaColony3456_20170921.m`
%             because realized that that code does not include the part to
%             rotate the modeled clicks.
%             The operation here is to generate modeled clicks at different
%             freq and noise levels, and will write another code to assemble
%             the composite modeled clicks
% 2017 10 03  modify the code to do multi-freq piston model simulation

clear
warning off

usrn = getenv('username');
if isunix
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/rbfinterp_v1.2');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/saveSameSize');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/EllipseDirectFit');
    data_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
    save_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
    model_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp/bp_bem_modeling/';
else
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    addpath('F:\Dropbox\0_CODE\MATLAB\EllipseDirectFit');
    data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    model_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp\bp_bem_modeling';
end

cvec = 0:-3:-39;

plot_opt_all_click = 1;
plot_opt_indiv_click = 1;
save_opt = 1;

results_path = 'analysis_results_figs';
indiv_data_path = 'rotate_all_click_2010308';
diff_data_path = 'model_bp_proj_RaColony_diffonly_20170921';
ss = strsplit(diff_data_path,'_');
model_shape = ss{end-2};
model_calc_date = ss{end};
indiv_data_file = dir(fullfile(data_base_path,results_path,indiv_data_path,'*.mat'));
diff_file = dir(fullfile(data_base_path,results_path,diff_data_path,'*.mat'));

A.param.model_base_path = model_base_path;
A.param.data_base_path = data_base_path;
A.param.diff_data_path = diff_data_path;
A.param.indiv_data_path = indiv_data_path;
A.param.indiv_data_file = indiv_data_file;


% Other params
freq_main = 35e3;  % freq used to make shift/rotation [Hz]
freq_other = [20:5:30,40:5:60]*1e3;  % other freq in simulation [Hz]
freq_all = [freq_main,freq_other];
noise_mean = 0;  % added noise profile [dB]
noise_std_all = [1,0,2];

A.freq.main = freq_main;
A.freq.other = freq_other;
A.freq.all = freq_all;
A.param.noise_mean = noise_mean;



%% Get important params from rotated data files
% Load 1 file to get rotation data info
D = load(fullfile(data_base_path,results_path,diff_data_path,diff_file(1).name));  % rotated data
num_ch = length(D.raw_meas_from_mic.az);
A.map = D.map;
A.param.rng_seed = D.param.rng_seed;
A.param.num_ch = num_ch;

% Set bounds for minvtran
[xxx,yyy] = mfwdtran(D.map.mstruct,[0,0,-90,90],[-180,180,0,0]);
xxx = xxx(1:2);
yyy = yyy(3:4);
xy_lim = [xxx(1) xxx(2) yyy(1) yyy(2)];



%% Load azel from experimental data to find best-fitting piston diameter
% Get info from individual files
for iF=1:length(indiv_data_file)
    % individual click data
    B_data = load(fullfile(data_base_path,results_path,...
                           indiv_data_path,indiv_data_file(iF).name));
    data(iF) = get_azel_params(B_data);
end

% Remove points out of xy map projection limit
data_idx_bad = [data(:).elps_x]<xy_lim(1) | [data(:).elps_x]>xy_lim(2) |...
               [data(:).elps_y]<xy_lim(3) | [data(:).elps_y]>xy_lim(4);

% Convert from x-y to az-el angles -- inverse map projection
[el_data,az_data] = minvtran(D.map.mstruct,[data(~data_idx_bad).elps_x]',...
                                           [data(~data_idx_bad).elps_y]');

% Set bp param
bp_info.c = 344;  % sound speed [m/s]
bp_info.freq = A.freq.main;  % [Hz]
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
match_deg = nanmean([az_data;el_data]);  % mean of az and el from exp data
[~,idx_mean] = min(abs(piston_3dB_c(:,1)-match_deg));
a_mean = piston_3dB_c(idx_mean,2);
close
bp_info.a = a_mean;
A.model_bp_info = bp_info;





for iN=1:length(noise_std_all)  % loop through all noise levels

noise_std = noise_std_all(iN);
A.param.noise_std = noise_std;

% Set up saving path
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(data_base_path,results_path,sprintf('%s_std%2.1f',script_name,noise_std));
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% re-seed the random number generator for each noise level
rng(D.param.rng_seed);

% Loop through all files
flag_trial = 0;
flag_count = 0;
for iS=1:length(diff_file)

    ss = strsplit(strtok(diff_file(iS).name,'.'),'_');
    ss2 = strsplit(script_name,'_');
    save_fname = strjoin([script_name,sprintf('std%2.1f',noise_std), ...
                        ss(end-2:end)],'_');
    title_str = sprintf('%s %dkHz noise std=%2.1f %s',...
                        ss2{4},freq_all(1)/1e3,noise_std,strjoin([ss(end-2:end)],'_'));
    
    % Load partial simulation results (bpctr diff)
    D = load(fullfile(data_base_path,results_path,diff_data_path, ...
                      diff_file(iS).name));
    disp(sprintf('file %03d: %s',iS,diff_file(iS).name));

    % Load rotate measurement file --> for checking/debugging
    %R = load(fullfile(data_base_path,results_path,D.param.rotate_data_path, ...
    %                  D.param.rotate_data_file));
    
    A.param.diff_data_file = diff_file(iS).name;

    % Set up figures
    if ~flag_trial
        t_name = strjoin(ss(7:8),'_');
        trial_file_all = dir(fullfile(data_base_path,...
                                      results_path,diff_data_path,['*_',t_name,'_*.mat']));
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
    
    for iF=1:length(freq_all)
        
        freq_model = freq_all(iF);
        disp(sprintf('Simulating freq = %dkHz',freq_model/1e3));

        % Model mic output

        % Simulation -- need to do for all freq
        if isempty(D.diff)  % if simulated data not available due to invalid
                            % minvtran values in raw measurements
            A.v_mic = [];
            A.click_side = [];
            A.raw = [];
            A.rot_max = [];
            A.rot_elpctr = [];
            A.rot_elpctr_tilt = [];
            
        else    % if simulated data available
                % Get beampattern model elq/azq
            v_mic = model_beam_dB(bp_info,...
                        [D.meas.ectr.az,D.meas.ectr.el]/180*pi,...  % all inputs in [rad]
                        [D.raw_meas_from_mic.az,D.raw_meas_from_mic.el],...
                        freq_model);

            % Add noise
            v_mic = v_mic+randn(size(v_mic))*noise_std+noise_mean;
            v_mic = v_mic-max(v_mic);

            A.v_mic(:,iF) = v_mic;  % save simulated mic recording


            % Only do shift/rotation at the first freq (35 kHz)
            if iF==1

                % Fit ellipse
                [click_side,raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
                    shift_rotate_bp(D.raw_meas_from_mic.az,D.raw_meas_from_mic.el,...
                                    v_mic,'eckert4',D.param.iterative_shift_threshold);

                % if the best-fitting ellipse for shift_max is outside of globe
                % then don't plot
                if isempty(rot_elpctr) || isempty(rot_elpctr_tilt)
                    flag_elps_fail = 1;
                else
                    flag_elps_fail = 0;
                end
                
                % Check az, el shift and tilt
                if plot_opt_indiv_click==1
                    Vpass.raw = raw;
                    Vpass.rot_max = rot_max;
                    Vpass.rot_elpctr = rot_elpctr;
                    Vpass.rot_elpctr_tilt = rot_elpctr_tilt;
                    
                    % Movement of [az,el] of each mic
                    fig_azel = plot_indiv_click_azel_movement(Vpass);
                    title(regexprep(title_str,'_','\\_'));
                    saveSameSize_res(fig_azel,120,'file',...
                                 fullfile(save_path,sprintf('%s_azel_chk.png',save_fname)),...
                                 'format','png','renderer','painters');
                    close(fig_azel)
                    
                    % Plot the shift/tilt procedure
                    fig_fit = plot_indiv_click_rotate(Vpass,cvec,D.map.mstruct);
                    suptitle(regexprep(title_str,'_','\\_'));
                    saveSameSize_res(fig_fit,120,'file',...
                                 fullfile(save_path,sprintf('%s_contour.png',save_fname)),...
                                 'format','png','renderer','painters');
                    close(fig_fit)
                end
                
                % Rotated model output
                A.click_side = click_side;
                A.raw = raw;
                A.rot_max = rot_max;
                A.rot_elpctr = rot_elpctr;
                A.rot_elpctr_tilt = rot_elpctr_tilt;
                
                % Plot rotated bp using ellipse center
                if flag_count<=length(trial_file_all) && flag_elps_fail~=1
                    if plot_opt_all_click
                        figure(fig_clicks)
                        subplot(numrow,4,flag_count);
                        axesm(D.map.mstruct)
                        contour(rot_elpctr_tilt.xq,rot_elpctr_tilt.yq,rot_elpctr_tilt.vq_norm,cvec,'fill','on');
                        framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
                        gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
                        tightmap
                        axis off
                        title(sprintf('Call #%02d',str2double(ss{end}(2:3))));
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
                        saveSameSize_res(fig_clicks,120,'file',...
                                     fullfile(save_path,sprintf('all_clicks_%s_%s.png',script_name,t_name)),...
                                     'format','png','renderer','painters');
                        close(fig_clicks)
                    end
                end
            end  % if at first freq
            
        end  % if values invalid for minvtran in raw measurement ellipse fitting

    end % loop through all freq

    % Plot to check all freq
    if ~isempty(A.raw)
        fig_all_freq_proj = figure('position',[100 100 900 650]);
        numrow_all_freq = ceil(length(freq_all)/3);
        for iF=1:length(freq_all)
            figure(fig_all_freq_proj);
            plot_bp_simple(subplot(numrow_all_freq,3,iF),A.raw.az,A.raw.el,A.v_mic(:,iF)', ...
                           A.map.map_projection);
            title(sprintf('%d kHz',freq_all(iF)/1e3));
        end
        figure(fig_all_freq_proj);
        suptitle(regexprep(title_str,'_','\\_'));
        saveSameSize_res(fig_all_freq_proj,120,'file',...
                         fullfile(save_path,sprintf('%s_all_freq_proj.png',save_fname)),...
                         'format','png','renderer','painters');
    end

    % Save output
    if save_opt==1
        save(fullfile(save_path,[save_fname,'.mat']),'-struct','A');
    end

end % loop through all diff_file

end % loop through all noise levels


warning on