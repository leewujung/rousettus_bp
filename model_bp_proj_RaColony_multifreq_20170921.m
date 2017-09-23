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
    data_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
    save_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
    model_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp/bp_bem_modeling/';
else
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
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

A.param.data_base_path = data_base_path;
A.param.diff_data_path = diff_data_path;

% Load 1 file to get info
D = load(fullfile(data_base_path,results_path,diff_data_path,diff_file(1).name));  % rotated data
num_ch = length(D.raw_meas_from_mic.az);
A.map = D.map;
A.param.num_ch = num_ch;

freq_model_all = [25:10:55]*1e3;
noise_mean = 0;  % added noise profile [dB]
noise_std = 1;

A.param.noise_mean = noise_mean;
A.param.noise_std = noise_std;
A.param.rng_seed = D.param.rng_seed;

% Set up saving path
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(data_base_path,results_path,sprintf('%s_std%2.1f',script_name,noise_std));
if ~exist(save_path,'dir')
    mkdir(save_path);
end


% Seed randomization
rng(D.param.rng_seed);  % seed the random number generator      

for iF=1:length(freq_model_all)
    
    freq_model = freq_model_all(iF);
    A.param.freq_model = freq_model;
    
    
    disp(sprintf('Simulating freq = %dkHz',freq_model/1e3));

    % Load simulated beampattern at freq_model
    sbp = strsplit(D.BP.bp_model_file,'_');
    bp_model_file = strjoin([sbp(1:6),sprintf('%dkHz',freq_model/1e3),sbp(8:end)],'_');
    BP = load(fullfile(model_base_path,D.BP.bp_model_path,bp_model_file));
    A.BP = BP;

    % Loop through all files
    flag_trial = 0;
    flag_count = 0;
    for iS=1:length(diff_file)

        ss = strsplit(strtok(diff_file(iS).name,'.'),'_');
        ss2 = strsplit(script_name,'_');
        save_fname = strjoin([script_name,sprintf('%dkHz_std%2.1f',freq_model/1e3,noise_std), ...
                            ss(end-2:end)],'_');
        title_str = strjoin([ss2(4:6),sprintf('%dkHz_std%2.1f',freq_model/ ...
                                              1e3,noise_std),ss(end-2:end)],'_');
        
        
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
        
        % Simulation
        if isempty(D.diff)  % if simulated data not available due to invalid
                            % minvtran values in raw measurements
            
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
            
            % Fit ellipse
            [click_side,raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
                shift_rotate_bp(D.raw_meas_from_mic.az,D.raw_meas_from_mic.el,...
                                v_mic,'eckert4',D.param.iterative_shift_threshold);
            
            % Check az, el shift and tilt
            if plot_opt_indiv_click==1
                Vpass.raw = raw;
                Vpass.rot_max = rot_max;
                Vpass.rot_elpctr = rot_elpctr;
                Vpass.rot_elpctr_tilt = rot_elpctr_tilt;
                
                % Movement of [az,el] of each mic
                fig_azel = plot_indiv_click_azel_movement(Vpass);
                title(regexprep(title_str,'_','\\_'));
                saveSameSize(fig_azel,'file',...
                             fullfile(save_path,sprintf('%s_azel_chk.png',save_fname)),...
                             'format','png','renderer','painters');
                close(fig_azel)
                
                % Plot the shift/tilt procedure
                fig_fit = plot_indiv_click_rotate(Vpass,cvec,D.map.mstruct);
                suptitle(regexprep(title_str,'_','\\_'));
                saveSameSize(fig_fit,'file',...
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
            if flag_count<=length(trial_file_all)
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
                    saveSameSize(fig_clicks,'file',...
                                 fullfile(save_path,sprintf('all_clicks_%s_%s.png',script_name,t_name)),...
                                 'format','png','renderer','painters');
                    close(fig_clicks)
                end
            end

            
        end  % if values invalid for minvtran in raw measurement ellipse fitting

        if save_opt==1
            save(fullfile(save_path,[save_fname,'.mat']),'-struct','A');
        end

        
    end % loop through all diff_file

end % loop through all freq