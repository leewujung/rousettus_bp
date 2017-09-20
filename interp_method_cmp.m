% 2016 07 22  Compare interpolation method

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

% noise_mean = 0;  % added noise profile [dB]
% noise_std = 0.5;   % [dB]
results_path = 'analysis_results_figs';
rotate_data_path = 'rotate_all_click_20160721';
match = 'azel';  % az|el|azel
res_deg = 10;     % interpolation resolution [deg]

% Find all click files to match
bat_proc_file_all = dir(fullfile(data_base_path,results_path,rotate_data_path,'*.mat'));

% Set various paths
[~,script_name,~] = fileparts(mfilename('fullpath'));
% script_name = sprintf('%s_%s_std%2.1f',script_name,match,noise_std);
save_path = fullfile(data_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

azel_distr_path = 'az_el_distr_model_data_cmp_20160722';
azel_distr_file = 'az_el_distr_model_data_cmp_20160722_data.mat';
azel_distr = load(fullfile(data_base_path,results_path,azel_distr_path,azel_distr_file));

S.base_path = data_base_path;
S.save_path = save_path;
S.azel_distr_path = azel_distr_path;
S.azel_distr_file = azel_distr_file;

% Param
freq_model = 35e3;
it_shift_th = 0.005;
azel_bnd = [300,150];  % bound of which the rotate-shift difference larger than this would be discarded
cvec = 0:-3:-39;

plot_opt_all_click = 1;
plot_opt_indiv_click = 1;
save_opt = 1;

S.freq_modeled = freq_model;
S.iterative_shift_threshold = it_shift_th;
S.shift_bound_azel = azel_bnd;


% Set map projection
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);

S.map_projection = map_proj;
S.mstruct = mstruct;


% Set bp param
bp_info.c = 344;  % sound speed [m/s]
bp_info.freq = freq_model;  % [Hz]
bp_info.k = 2*pi*bp_info.freq/bp_info.c;  % wavenumber
bp_info.type = 'piston';


% Find aperture with best-fitting -3dB contour to mean azimuth of data
a_fine = (5:0.01:7)*1e-3; % [m]
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


% Simualte for each click
for iB = 1:length(bat_proc_file_all)
    bat_proc_file = bat_proc_file_all(iB).name;
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

    [~,rbf_vq_norm,azq,elq] = interp_bp_res(az,el,call_dB,'rbf',res_deg);
    [~,nn_vq_norm,~,~] = interp_bp_res(az,el,call_dB,'natural',res_deg);
    
    gt_bp = model_beam_dB(bp_info,call_max_azel,[azq(:) elq(:)]);  % ground-truth bp
    gt_bp = reshape(gt_bp,size(azq));
    
    az = az/pi*180;  % convert to [deg]
    el = el/pi*180;
    azq = azq/pi*180;
    elq = elq/pi*180;
    
    % Find boundary
    k = boundary(az,el,0);  % outer boundary of all measured points
    [in,on] = inpolygon(azq,elq,az(k),el(k));
    in = in|on;
    
    % Set values outside of az-el boundary to NaN
    gt_bp(~in) = NaN;
    rbf_vq_norm(~in) = NaN;
    nn_vq_norm(~in) = NaN;

    % Save errors
    nn_err{iB} = reshape(gt_bp-nn_vq_norm,[],1);
    rbf_err{iB} = reshape(gt_bp-rbf_vq_norm,[],1);
end

nn_err_all = cell2mat(nn_err');
rbf_err_all = cell2mat(rbf_err');

mx = max([rbf_err_all;nn_err_all]);
mn = min([rbf_err_all;nn_err_all]);

figure
subplot(211)
hist(rbf_err_all,mn:0.1:mx);
legend('Radial basis fcn')
grid
subplot(212)
hist(nn_err_all,mn:0.1:mx);
legend('Natural beighbor')
grid

nanmean([rbf_err_all, nn_err_all])
nanmedian([rbf_err_all, nn_err_all])


