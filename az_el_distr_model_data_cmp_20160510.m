% 2016 04 19  Plot distribution of azimuth and elevation for best-fitting
%             ellipse
% 2016 05 10  Revise and plot for paper

clear
warning off
usrn = getenv('username');

if strcmp(usrn,'Wu-Jung')
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    addpath('F:\Dropbox\0_CODE\MATLAB\swtest');
else
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\swtest']);
end

% Set up various paths
[~,script_name,~] = fileparts(mfilename('fullpath'));
if strcmp(usrn,'Wu-Jung')
    base_path = ['F:\Dropbox\0_ANALYSIS\bp_processing'];
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
end

save_root = 'analysis_results_figs';
save_path = fullfile(base_path,save_root,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

data_path = 'rotate_all_click_20160419';
data_file_all = dir(fullfile(base_path,save_root,data_path,'*.mat'));
model_path = 'model_piston_proj_20160510';
model_file_all = dir(fullfile(base_path,save_root,model_path,'*.mat'));

A.data_path = data_path;
A.data_file_all = data_file_all;
A.model_path = model_path;
A.model_file_all = model_file_all;

data_file = 'rotate_all_click_20160419_all_click_rotated_3bat.mat';
model_file = 'model_piston_proj_20160510_all_model_rotated_3bats.mat';

A.data_file = data_file;
A.model_file = model_file;

% for iF=1%:length(data_file_all)
    B_data = load(fullfile(base_path,save_root,data_path,data_file));
    B_model = load(fullfile(base_path,save_root,model_path,model_file));
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
    theta = nan(trial_num*ccmm,2);
    for iT=1:trial_num
        for iC=1:click_num(iT)
            E_data = B_data.rotate_data{iT}(iC).rot_elpctr_tilt.E;
            E_model = B_model.rotate_data{iT}(iC).rot_elpctr_tilt.E;
            a0((iT-1)*ccmm+iC,1) = E_data.a0;
            b0((iT-1)*ccmm+iC,1) = E_data.b0;
            e((iT-1)*ccmm+iC,1) = E_data.e;
            ar((iT-1)*ccmm+iC,1) = E_data.ar;
            theta((iT-1)*ccmm+iC,1) = E_data.theta;
            a0((iT-1)*ccmm+iC,2) = E_model.a0;
            b0((iT-1)*ccmm+iC,2) = E_model.b0;
            e((iT-1)*ccmm+iC,2) = E_model.e;
            ar((iT-1)*ccmm+iC,2) = E_model.ar;
            theta((iT-1)*ccmm+iC,2) = E_model.theta;
        end
    end
% end
nanidx = isnan(a0(:,1));
a0(nanidx,:) = [];
b0(nanidx,:) = [];
e(nanidx,:) = [];
ar(nanidx,:) = [];
theta(nanidx,:) = [];

% Fix major/minor axes
idx_data = a0(:,1)>b0(:,1);
minor_data = a0(:,1);
major_data = b0(:,1);
minor_data(idx_data) = b0(idx_data,1);
major_data(idx_data) = a0(idx_data,1);

A.minor_data = minor_data;
A.major_data = major_data;

idx_model = a0(:,2)>b0(:,2);
minor_model = a0(:,2);
major_model = b0(:,2);
minor_model(idx_model) = b0(idx_model,2);
major_model(idx_model) = a0(idx_model,2);

A.minor_model = minor_model;
A.major_model = major_model;

% Convert from x-y to az-el angles
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);
[el_data,az_data] = minvtran(mstruct,minor_data,major_data);
[el_model,az_model] = minvtran(mstruct,minor_model,major_model);

A.az = az_data;
A.el = el_data;
A.az = az_model;
A.el = el_model;

% Output mean and std
az_mean_data = mean(az_data);
el_mean_data = mean(el_data);
az_std_data = std(az_data);
el_std_data = std(el_data);
az_mean_model = mean(az_model);
el_mean_model = mean(el_model);
az_std_model = std(az_model);
el_std_model = std(el_model);

A.az_mean_data = az_mean_data;
A.el_mean_data = el_mean_data;
A.az_std_data = az_std_data;
A.el_std_data = el_std_data;
A.az_mean_model = az_mean_model;
A.el_mean_model = el_mean_model;
A.az_std_model = az_std_model;
A.el_std_model = el_std_model;

save(fullfile(save_path,[script_name,'_azel_xy_distr.mat']),'-struct','A');


% Plot: summary of data
fig = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(231)
plot(minor_data,major_data,'.')
xlabel('X projection')
ylabel('Y projection')
title('X-Y')
axis equal
grid
axis([0 1.5 0 1.5])

subplot(232)
hist(minor_data,0:0.05:1.5)
xlabel('X projection');
ylabel('Count')
title('X distr')

subplot(233)
hist(major_data,0:0.05:1.5)
xlabel('Y projection');
ylabel('Count')
title('Y distr')

subplot(234)
plot(az_data,el_data,'.')
xlabel('Azimuth (deg)')
ylabel('Elevation (deg)')
title('AZ-EL')
axis equal
grid
axis([0 90 0 90])

subplot(235)
hist(az_data,0:5:90)
xlabel('Azimuth (deg)')
ylabel('Count')
title('AZ distr')

subplot(236)
hist(el_data,0:5:90)
xlabel('Elevation (deg)')
ylabel('Count')
title('EL distr')

suptitle('Data summary');

saveas(fig,fullfile(save_path,[script_name,'_azel_xy_distr_data.fig']),'fig');
saveSameSize(fig,'file',...
    fullfile(save_path,[script_name,'_azel_xy_distr_data']),...
    'format','png','renderer','painters');

% Plot: summary of model output
fig = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(231)
plot(minor_model,major_model,'.')
xlabel('X projection')
ylabel('Y projection')
title('X-Y')
axis equal
grid
axis([0 1.5 0 1.5])

subplot(232)
hist(minor_model,0:0.05:1.5)
xlabel('X projection');
ylabel('Count')
title('X distr')

subplot(233)
hist(major_model,0:0.05:1.5)
xlabel('Y projection');
ylabel('Count')
title('Y distr')

subplot(234)
plot(az_model,el_model,'.')
xlabel('Azimuth (deg)')
ylabel('Elevation (deg)')
title('AZ-EL')
axis equal
grid
axis([0 90 0 90])

subplot(235)
hist(az_model,0:5:90)
xlabel('Azimuth (deg)')
ylabel('Count')
title('AZ distr')

subplot(236)
hist(el_model,0:5:90)
xlabel('Elevation (deg)')
ylabel('Count')
title('EL distr')

suptitle('Model summary');

saveas(fig,fullfile(save_path,[script_name,'_azel_xy_distr_model.fig']),'fig');
saveSameSize(fig,'file',...
    fullfile(save_path,[script_name,'_azel_xy_distr_model']),...
    'format','png','renderer','painters');

% Plot: data-model comparison
fig = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(231)
plot(minor_data,major_data,'.',minor_model,major_model,'.')
xlabel('X projection')
ylabel('Y projection')
title('X-Y')
axis equal
grid
axis([0 1.5 0 1.5])
l = legend('Data','Model');
set(l,'fontsize',12)

subplot(232)
hist([minor_data,minor_model],0:0.05:1.5)
xlabel('X projection');
ylabel('Count')
title('X distr')
l = legend('Data','Model');
set(l,'fontsize',12)

subplot(233)
hist([major_data,major_model],0:0.05:1.5)
xlabel('Y projection');
ylabel('Count')
title('Y distr')

subplot(234)
plot(az_data,el_data,'.',az_model,el_model,'.')
xlabel('Azimuth (deg)')
ylabel('Elevation (deg)')
title('AZ-EL')
axis equal
grid
axis([0 90 0 90])

subplot(235)
hist([az_data,az_model],0:5:90)
xlabel('Azimuth (deg)')
ylabel('Count')
title('AZ distr')

subplot(236)
hist([el_data,el_model],0:5:90)
xlabel('Elevation (deg)')
ylabel('Count')
title('EL distr')

suptitle('Data-model comparison');

saveas(fig,fullfile(save_path,[script_name,'_azel_xy_distr_data_model_cmp.fig']),'fig');
saveSameSize(fig,'file',...
    fullfile(save_path,[script_name,'_azel_xy_distr_data_model_cmp']),...
    'format','png','renderer','painters');

% Plot: el/az ratio data-model comparison
fig = figure;
nedges = 0:0.2:4;
nedges_center = nedges+diff(nedges(1:2))/2;
bc_data = histc(el_data./az_data,nedges);
bc_model = histc(el_model./az_model,nedges);

subplot(211)
bar(nedges_center,bc_data/sum(bc_data));
ylabel('Relative frequency');
l = legend('Data');
set(l,'fontsize',12)

subplot(212)
bar(nedges_center,bc_model/sum(bc_model));
ylabel('Relative frequency');
xlabel('az/el ratio');
l = legend('Model');
set(l,'fontsize',12)

saveas(fig,fullfile(save_path,[script_name,'_azel_ratio_distr_data_model_cmp.fig']),'fig');
saveSameSize(fig,'file',...
    fullfile(save_path,[script_name,'_azel_ratio_distr_data_model_cmp']),...
    'format','png','renderer','painters');

[H_data,p_data,~] = swtest(el_data./az_data,0.1);
[H_model,p_model,~] = swtest(el_model./az_model,0.1);

