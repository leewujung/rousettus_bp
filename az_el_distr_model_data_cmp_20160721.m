% 2016 04 19  Plot distribution of azimuth and elevation for best-fitting
%             ellipse
% 2016 05 10  Revise and plot for paper

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

data_path = 'rotate_all_click_20160419';
model_path = 'model_piston_proj_20160510';

A.data_path = data_path;
A.model_path = model_path;

data_file = 'rotate_all_click_20160419_all_click_rotated_3bat.mat';
model_file = 'model_piston_proj_20160510_all_model_rotated_3bats.mat';

A.data_file = data_file;
A.model_file = model_file;


% Load data/model
B_data = load(fullfile(data_base_path,results_path,data_path,data_file));
B_model = load(fullfile(data_base_path,results_path,model_path,model_file));
trial_num = length(B_data.rotate_data);

click_num = zeros(trial_num,1);
for iT=1:trial_num
    click_num(iT) = length(B_data.rotate_data{iT});
end
ccmm = max(click_num);
a0 = nan(trial_num*ccmm,2);
b0 = nan(trial_num*ccmm,2);
e = nan(trial_num*ccmm,2);
ar = nan(trial_num*ccmm,2);
elps_x_max = nan(trial_num*ccmm,2);
elps_y_max = nan(trial_num*ccmm,2);
theta = nan(trial_num*ccmm,2);
for iT=1:trial_num
    for iC=1:click_num(iT)
        E_data = B_data.rotate_data{iT}(iC).rot_elpctr_tilt.E;
        [ar_tmp,elps_xy_tmp] = get_ellipse_ar(E_data);
        ar((iT-1)*ccmm+iC,1) = ar_tmp;
        elps_x_max((iT-1)*ccmm+iC,1) = range(elps_xy_tmp(:,1))/2;
        elps_y_max((iT-1)*ccmm+iC,1) = range(elps_xy_tmp(:,2))/2;

        E_model = B_model.rotate_data{iT}(iC).rot_elpctr_tilt.E;
        [ar_tmp,elps_xy_tmp] = get_ellipse_ar(E_model);
        ar((iT-1)*ccmm+iC,2) = ar_tmp;
        elps_x_max((iT-1)*ccmm+iC,2) = range(elps_xy_tmp(:,1))/2;
        elps_y_max((iT-1)*ccmm+iC,2) = range(elps_xy_tmp(:,2))/2;
        
        % NOTE: ar is calculated again here separately because the one
        % comes in E assumes whichever longer is along the Y axis. This
        % may not be true for piston model simulation results.
        
        a0((iT-1)*ccmm+iC,1) = E_data.a0;
        b0((iT-1)*ccmm+iC,1) = E_data.b0;
        e((iT-1)*ccmm+iC,1) = E_data.e;
        %ar((iT-1)*ccmm+iC,1) = E_data.ar;
        theta((iT-1)*ccmm+iC,1) = E_data.theta;
        
        a0((iT-1)*ccmm+iC,2) = E_model.a0;
        b0((iT-1)*ccmm+iC,2) = E_model.b0;
        e((iT-1)*ccmm+iC,2) = E_model.e;
        %ar((iT-1)*ccmm+iC,2) = E_model.ar;
        theta((iT-1)*ccmm+iC,2) = E_model.theta;
    end
end
nanidx = isnan(a0(:,1));
a0(nanidx,:) = [];
b0(nanidx,:) = [];
e(nanidx,:) = [];
ar(nanidx,:) = [];
theta(nanidx,:) = [];
elps_x_max(nanidx,:) = [];
elps_y_max(nanidx,:) = [];

A.a0 = a0;
A.b0 = b0;
A.e = e;
A.ar = ar;
A.theta = theta;
A.elps_x_max = elps_x_max;
A.elps_y_max = elps_y_max;


% Convert from x-y to az-el angles
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);
[el_data,az_data] = minvtran(mstruct,elps_x_max(:,1),elps_y_max(:,1));
[el_model,az_model] = minvtran(mstruct,elps_x_max(:,2),elps_y_max(:,2));

A.el_data = el_data;
A.az_data = az_data;
A.el_model = el_model;
A.az_model = az_model;

save(fullfile(save_path,[script_name,'_azel_xy.mat']),'-struct','A');

% Plot
fig_az_el = figure;
plot(az_data,el_data,'.');
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
hist(el_data./az_data,0:0.1:4);
subplot(212)
hist(el_model./az_model,0:0.1:4);

save_fname = sprintf('%s_ar',script_name);
saveas(fig_ar,...
    fullfile(save_path,[save_fname,'.fig']),'fig');
saveSameSize(fig_ar,'file',fullfile(save_path,[save_fname,'.png']),...
    'format','png','renderer','painters');
