% 2016 04 19  Plot distribution of azimuth and elevation for best-fitting
%             ellipse
% 2016 05 10  Revise and plot for paper
% 2015 07 22  Update routines to load data and model parts
% 2015 07 23  Load models with different SNR

clear
warning off

if isunix
    addpath('~/Dropbox/0_CODE/MATLAB/saveSameSize');
    addpath('~/Dropbox/0_CODE/MATLAB/swtest');

    data_base_path = '~/internal_2tb/Dropbox/0_ANALYSIS/bp_processing';
    save_base_path = '~/Dropbox/Z_wjlee/projects/rousettus_bp';
else
    usrn = getenv('username');
    if strcmp(usrn,'Wu-Jung')
        addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
        addpath('F:\Dropbox\0_CODE\MATLAB\swtest');
    else
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\swtest']);
    end

    if strcmp(usrn,'Wu-Jung')
        data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
        save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    else
        data_base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
        save_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
    end
end


% Set up various paths
results_path = 'analysis_results_figs';

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

model_path = 'model_piston_proj_20160721_azel_std0.5';
model_file = dir(fullfile(data_base_path,results_path,model_path,'*.mat'));
A.model_path = model_path;
A.model_file = model_file;

for iF=1:length(model_file)
    B_model = load(fullfile(data_base_path,results_path,...
            model_path,model_file(iF).name));
    [ar_tmp,elps_xy_tmp] = get_ellipse_ar(B_model.rot_elpctr_tilt.E);
    ar_model(iF,1) = ar_tmp;
    elps_x_max_model(iF,1) = range(elps_xy_tmp(:,1))/2;
    elps_y_max_model(iF,1) = range(elps_xy_tmp(:,2))/2;
    e_model(iF,:) = B_model.rot_elpctr_tilt.E.e;
    a0_model(iF,:) = B_model.rot_elpctr_tilt.E.a0;
    b0_model(iF,:) = B_model.rot_elpctr_tilt.E.b0;
    theta_model(iF,:) = B_model.rot_elpctr_tilt.E.theta;
end

A.model.ar = ar_model;
A.model.a0 = a0_model;
A.model.b0 = b0_model;
A.model.e = e_model;
A.model.theta = theta_model;
A.model.elps_x_max = elps_x_max_model;
A.model.elps_y_max = elps_y_max_model;


% Convert from x-y to az-el angles
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);
[el_model,az_model] = minvtran(mstruct,elps_x_max_model(:,1),elps_y_max_model(:,1));
% [el_model,az_model] = minvtran(mstruct,elps_x_max(:,2),elps_y_max(:,2));

A.el_data = el_model;
A.az_data = az_model;
% A.el_model = el_model;
% A.az_model = az_model;

save(fullfile(save_path,[script_name,'_model.mat']),'-struct','A');

% Plot
fig_az_el = figure;
plot(az_model,el_model,'.');
hold on
plot(az_model,el_model,'x');
xlabel('Azimuth');
ylabel('Elevation');
legend('Data','Model');
grid
axis equal

save_fname = sprintf('%s_az_el_scatter',script_name);
saveas(fig_az_el,...
    fullfile(save_path,[save_fname,'.fig']),'fig');
saveSameSize(fig_az_el,'file',fullfile(save_path,[save_fname,'.png']),...
    'format','png','renderer','painters');


fig_ar = figure;
subplot(211)
hist(el_model./az_model,0:0.1:4);
subplot(212)
hist(el_model./az_model,0:0.1:4);

save_fname = sprintf('%s_ar',script_name);
saveas(fig_ar,...
    fullfile(save_path,[save_fname,'.fig']),'fig');
saveSameSize(fig_ar,'file',fullfile(save_path,[save_fname,'.png']),...
    'format','png','renderer','painters');
