% 2016 04 19  Plot distribution of azimuth and elevation for best-fitting
%             ellipse

clear
warning off
usrn = getenv('username');

if strcmp(usrn,'Wu-Jung')
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    addpath('F:\Dropbox\0_CODE\MATLAB\export_fig-master');
else
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\export_fig-master']);
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

A.data_path = data_path;
A.data_file_all = data_file_all;

for iF=1:length(data_file_all)
    B = load(fullfile(base_path,save_root,data_path,data_file_all(1).name));
    trial_num = length(B.rotate_data);
    
    for iT=1:trial_num
        click_num(iT) = length(B.rotate_data{iT});
    end
    ccmm = max(click_num);
    a0 = nan(trial_num*ccmm,1);
    b0 = nan(trial_num*ccmm,1);
    e = nan(trial_num*ccmm,1);
    ar = nan(trial_num*ccmm,1);
    ar_a0_b0 = nan(trial_num*ccmm,1);
    theta = nan(trial_num*ccmm,1);
    for iT=1:trial_num
        for iC=1:click_num(iT)
            E = B.rotate_data{iT}(iC).rot_elpctr_tilt.E;
            a0((iT-1)*ccmm+iC) = E.a0;
            b0((iT-1)*ccmm+iC) = E.b0;
            e((iT-1)*ccmm+iC) = E.e;
            ar((iT-1)*ccmm+iC) = E.ar;
            ar_a0_b0((iT-1)*ccmm+iC) = E.a0/E.b0;
            theta((iT-1)*ccmm+iC) = E.theta;
        end
    end
end
nanidx = isnan(a0);
a0(nanidx) = [];
b0(nanidx) = [];
e(nanidx) = [];
ar(nanidx) = [];
ar_a0_b0(nanidx) = [];
theta(nanidx) = [];

% Fix major/minor axes
idx = a0>b0;
minor = a0;
major = b0;
minor(idx) = b0(idx);
major(idx) = a0(idx);

A.minor = minor;
A.major = major;

% Convert from x-y to az-el angles
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);
[el,az] = minvtran(mstruct,minor,major);

A.az = az;
A.el = el;

% Output mean and std
az_mean = mean(az);
el_mean = mean(el);
az_std = std(az);
el_std = std(el);

A.az_mean = az_mean;
A.el_mean = el_mean;
A.az_std = az_std;
A.el_std = el_std;

save(fullfile(save_path,[script_name,'_azel_xy_distr.mat']),'-struct','A');


% Plot
fig = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(231)
plot(minor,major,'o')
xlabel('X projection')
ylabel('Y projection')
title('X-Y')
axis equal
grid
axis([0 1.5 0 1.5])

subplot(232)
hist(minor,0:0.05:1.5)
xlabel('X projection');
ylabel('Count')
title('X distr')

subplot(233)
hist(major,0:0.05:1.5)
xlabel('Y projection');
ylabel('Count')
title('Y distr')

subplot(234)
plot(az,el,'o')
xlabel('Azimuth (deg)')
ylabel('Elevation (deg)')
title('AZ-EL')
axis equal
grid
axis([0 90 0 90])

subplot(235)
hist(az,0:5:90)
xlabel('Azimuth (deg)')
ylabel('Count')
title('AZ distr')

subplot(236)
hist(el,0:5:90)
xlabel('Elevation (deg)')
ylabel('Count')
title('EL distr')

saveSameSize(fig,'file',...
    fullfile(save_path,[script_name,'_azel_xy_distr']),...
    'format','png','renderer','painters');

