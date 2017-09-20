% 2016 07 21  Try to plot the best-fitting ellipse in different map
%             projections 


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
    % Set up various paths
    if strcmp(usrn,'Wu-Jung')
        data_base_path = 'F:\Dropbox\0_ANALYSIS\bp_processing';
        save_base_path = 'F:\Dropbox\0_ANALYSIS\bp_processing';
    else
        data_base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
        save_base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
    end
end

% Make new folder for results from this code
[~,script_name,~] = fileparts(mfilename('fullpath'));
results_path = 'analysis_results_figs';
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

data_path = 'rotate_all_click_20160419';
data_file_all = dir(fullfile(data_base_path,results_path,data_path,'*.mat'));

A.data_path = data_path;
A.data_file_all = data_file_all;

data_file = 'rotate_all_click_20160419_all_click_rotated_39184.mat';

A.data_file = data_file;

% Load data
B_data = load(fullfile(data_base_path,results_path,data_path,data_file));

iT = 1;
iC = 1;
E = B_data.rotate_data{iT}(iC).rot_elpctr_tilt.E;

c3db_xy = E.c3db_xy;
xmin = min(c3db_xy(:,1))-range(c3db_xy(:,1)*0.5);
xmax = max(c3db_xy(:,1))+range(c3db_xy(:,1)*0.5);
ymin = min(c3db_xy(:,2))-range(c3db_xy(:,1)*0.5);
ymax = max(c3db_xy(:,2))+range(c3db_xy(:,1)*0.5);

fit_df = ezplot(E.eqt,[xmin,xmax,ymin,ymax]);  % plot ellipse
tmp = get(fit_df,'contourMatrix');

% best-fittign ellipse in eckert 4 projection
elps_x = tmp(1,2:end)';
elps_y = tmp(2,2:end)';  % best-fitting ellipse y coord

figure;
axesm('eckert4')
plot(elps_x,elps_y,'linewidth',2);
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
tightmap
title('Plot using (x,y) from rotation');

% transform to map (el,az) coordinate
axesm('eckert4')
[elps_el,elps_az] = minvtran(elps_x,elps_y);  

figure;
axesm('eckert4')
plotm(elps_el,elps_az,'linewidth',2);
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
tightmap
title('Plot using (el,az) after minvtran');


% transform to orthographic projection
axesm('ortho')
[elps_x_ortho,elps_y_ortho] = mfwdtran(elps_el,elps_az);

figure;
axesm('ortho')
plot(elps_x_ortho,elps_y_ortho,'linewidth',2);
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
tightmap
title('Plot using (x,y) when inv/fwd tran re ortho');


