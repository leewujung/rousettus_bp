% 2016 08 08  Assemble composite clicks from simulated mic receptions
% 2016 10 25  Update for version 1025 with out-of-bound points
% 2017 09 22  Modififed from `model_composite_RaColony3456_20170921.m`
%             because realized that that code does not include the part to
%             rotate the modeled clicks.
%             The operation here is to generate modeled clicks at different
%             freq and noise levels, and will write another code to assemble
%             the composite modeled clicks

clear

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
diff_data_path = 'model_bp_proj_RaColony_diffonly_20170921';
ss = strsplit(diff_data_path,'_');
model_shape = ss{end-2};
model_calc_date = ss{end};
diff_file = dir(fullfile(data_base_path,results_path,diff_data_path,'*.mat'));

A.param.model_base_path = model_base_path;
A.param.data_base_path = data_base_path;
A.param.diff_data_path = diff_data_path;

% Load 1 file to get info
D = load(fullfile(data_base_path,results_path,diff_data_path,diff_file(1).name));  % rotated data
num_ch = length(D.raw_meas_from_mic.az);
A.map = D.map;
A.param.num_ch = num_ch;

freq_main = 35e3;  % freq used to make shift/rotation [Hz]
freq_other = [20:5:30,40:5:60]*1e3;  % other freq in simulation [Hz]
freq_all = [freq_main,freq_other];
noise_mean = 0;  % added noise profile [dB]
noise_std_all = [1,0,2];

A.freq.main = freq_main;
A.freq.other = freq_other;
A.freq.all = freq_all;
A.param.noise_mean = noise_mean;
A.param.rng_seed = D.param.rng_seed;


for iN=1:length(noise_std_all)

noise_std = noise_std_all(iN);
A.param.noise_std = noise_std;

% Set up saving path
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(data_base_path,results_path,sprintf('%s_std%2.1f',script_name,noise_std));
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Seed randomization
rng(D.param.rng_seed);  % seed the random number generator      

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
    
    % Load partial simulation results
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

        % Load simulated beampattern at freq_model
        sbp = strsplit(D.BP.bp_model_file,'_');
        bp_model_file = strjoin([sbp(1:6),sprintf('%dkHz',freq_model/1e3),sbp(8:end)],'_');
        BP = load(fullfile(model_base_path,D.BP.bp_model_path,bp_model_file));
        BP.bp_model_path = D.BP.bp_model_path;
        BP.bp_model_file = bp_model_file;
        A.BP(iF) = BP;

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
            if D.raw_meas_from_mic.click_side==1  % right click --> no need to flip az
                az_model = BP.az;
            else   % left click --> need to flip az
                az_model = -BP.az;
            end
            el_model = BP.el;
            
            % Rotate model bp according to el_diff/az_diff from 35 kHz
            [el_model_rot,az_model_rot] = rotatem(el_model,az_model,...
                                                  [D.diff.el,D.diff.az],...
                                                  'inverse','degrees');

            % Project model bp to mic loc
            idxnotnan = ~isnan(BP.pp_plot);
            v_mic = rbfinterp([D.raw_meas_from_mic.az(:)'/pi*180;D.raw_meas_from_mic.el(:)'/pi*180],...
                              rbfcreate([az_model_rot(idxnotnan)';el_model_rot(idxnotnan)'],...
                                        BP.pp_plot(idxnotnan)',...
                                        'RBFFunction','multiquadrics'));

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
        fig_all_freq_proj = figure('position',[140 100 600 480]);
        fig_all_freq_bp = figure('position',[140 100 600 480]);;
        numrow_all_freq = ceil(length(freq_all)/2);
        for iF=1:length(freq_all)
            figure(fig_all_freq_proj);
            plot_bp_simple(subplot(numrow_all_freq,2,iF),A.raw.az,A.raw.el,A.v_mic(:,iF)', ...
                           A.map.map_projection);
            title(sprintf('%d kHz',freq_all(iF)/1e3));
            figure(fig_all_freq_bp);
            plot_bp_simple(subplot(numrow_all_freq,2,iF),A.BP(iF).az,A.BP(iF).el, ...
                           A.BP(iF).pp_plot,A.map.map_projection);
            title(sprintf('%d kHz',freq_all(iF)/1e3));
        end
        figure(fig_all_freq_proj);
        suptitle(regexprep(title_str,'_','\\_'));
        saveSameSize_res(fig_all_freq_proj,120,'file',...
                         fullfile(save_path,sprintf('%s_all_freq_proj.png',save_fname)),...
                         'format','png','renderer','painters');
        close(fig_all_freq_proj)
        figure(fig_all_freq_bp);
        suptitle(regexprep(title_str,'_','\\_'));
        saveSameSize_res(fig_all_freq_bp,120,'file',...
                         fullfile(save_path,sprintf('%s_all_freq_bp.png',save_fname)),...
                         'format','png','renderer','painters');
        close(fig_all_freq_bp)
    end

    % Save output
    if save_opt==1
        save(fullfile(save_path,[save_fname,'.mat']),'-struct','A');
    end

end % loop through all diff_file

end % loop through all noise levels