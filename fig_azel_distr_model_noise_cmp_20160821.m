% 2016 08 09  Plot distribution of azimuth and elevation for best-fitting
%             ellipse, for individual clicks and composite clicks

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
        data_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
        save_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
    end
end

results_path = 'analysis_results_figs';
model_noise_level = 0:2;

% Set up various paths
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

for iN=1:length(model_noise_level)
    
    % Path to model
    model_path = sprintf('model_piston_proj_20160818_azel_std%2.1f',model_noise_level(iN));
    model_file = dir(fullfile(data_base_path,results_path,model_path,'*.mat'));
    save_fname_iN = sprintf('%s_nstd%2.1f_data.mat',script_name,model_noise_level(iN));
    
    if exist(fullfile(save_path,save_fname_iN),'file')
        continue;
    end
    
    A.model_noise_level = model_noise_level(iN);
    A.model_path = model_path;
    A.model_file = model_file;
    
    % Loop through all files
    for iF=1:length(model_file)
        disp(['Loading ',model_file(iF).name]);
        B_model = load(fullfile(data_base_path,results_path,...
            model_path,model_file(iF).name));
        if ~isempty(B_model.rot_elpctr_tilt)
            [ar_tmp,elps_xy_tmp] = get_ellipse_ar(B_model.rot_elpctr_tilt.E);
            ar_model(iF,1) = ar_tmp;
            elps_x_max_model(iF,1) = range(elps_xy_tmp(:,1))/2;
            elps_y_max_model(iF,1) = range(elps_xy_tmp(:,2))/2;
            e_model(iF,:) = B_model.rot_elpctr_tilt.E.e;
            a0_model(iF,:) = B_model.rot_elpctr_tilt.E.a0;
            b0_model(iF,:) = B_model.rot_elpctr_tilt.E.b0;
            theta_model(iF,:) = B_model.rot_elpctr_tilt.E.theta;
        else
            ar_model(iF,1) = NaN;
            elps_x_max_model(iF,1) = NaN;
            elps_y_max_model(iF,1) = NaN;
            e_model(iF,:) = NaN;
            a0_model(iF,:) = NaN;
            b0_model(iF,:) = NaN;
            theta_model(iF,:) = NaN;
        end
    end
    A.ar_model = ar_model;
    A.elps_x_max_model = elps_x_max_model;
    A.elps_y_max_model = elps_y_max_model;
    A.e_model = e_model;
    A.a0_model = a0_model;
    A.b0_model = b0_model;
    A.theta_model = theta_model;
    
    % Convert from x-y to az-el angles
    map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
    mstruct = defaultm(map_proj);
    mstruct = defaultm(mstruct);
    
    % Revmoe points out of xy map projection limit
    [xxx,yyy] = mfwdtran(mstruct,[0,0,-90,90],[-180,180,0,0]);
    xxx = xxx(1:2);
    yyy = yyy(3:4);
    xy_lim = [xxx(1) xxx(2) yyy(1) yyy(2)];
    idx_bad_model = elps_x_max_model<xy_lim(1) | elps_x_max_model>xy_lim(2) |...
        elps_y_max_model<xy_lim(3) | elps_y_max_model>xy_lim(4);
    
    A.idx_bad_model = idx_bad_model;
    
    % Inverse map projection
    [el_model,az_model] = minvtran(mstruct,...
        elps_x_max_model(~idx_bad_model),elps_y_max_model(~idx_bad_model));
    
    A.el_model = el_model;
    A.az_model = az_model;
    
    % Save data
    save(fullfile(save_path,save_fname_iN),'-struct','A');
    
end % Model noise level


% For whatever reason the figures got closed when loading new file...
% so use the following to circumvent the problem
for iN=length(model_noise_level):-1:1
    save_fname_iN = sprintf('%s_nstd%2.1f_data.mat',script_name,model_noise_level(iN));
    A{iN} = load(fullfile(save_path,save_fname_iN));
end

fig_az_el = figure('position',[20 250 1300 420]);
corder = get(gca,'colororder');
cgrey = 200/255*[1 1 1];
for iN=1:length(model_noise_level)
    subplot(1,length(model_noise_level),iN)
    plot([0 90],[0 90],'color',cgrey);
    hold on
    plot(A{iN}.az_model(~A{iN}.idx_bad_model),...
        A{iN}.el_model(~A{iN}.idx_bad_model),'.',...
        'color',corder(4-iN,:),'markersize',6);
    axis equal
    axis([0 60 0 60])
    xlabel('Azimuth','fontsize',14)
    ylabel('Elevation','fontsize',14)
    set(gca,'fontsize',12,'xtick',0:10:90,'ytick',0:10:90)
    grid
    title(sprintf('Noise std %2.1f',model_noise_level(iN)));
end

saveas(fig_az_el,...
    fullfile(save_path,[script_name,'_scatter.fig']),'fig');
saveSameSize(fig_az_el,'file',fullfile(save_path,[script_name,'_scatter.png']),...
    'format','png','renderer','painters');
epswrite(fullfile(save_path,[script_name,'_scatter.eps']));


% % Plot aspect ratio histogram
% fig_ar = figure;
% subplot(211)
% hd = histogram(el_data(~idx_bad)./az_data(~idx_bad),...
%     0:0.1:4,'normalization','probability');
% hold on
% hleft = plot([1 1]*el_comp_left/az_comp_left,[0 1]);
% hright = plot([1 1]*el_comp_right/az_comp_right,[0 1]);
% ylim([0 0.18])
% title(title_text)
% subplot(212)
% hm = histogram(el_model(~idx_bad_model)./az_model(~idx_bad_model),...
%     0:0.1:4,'normalization','probability');
% xlabel('Ratio el/az');
% ylabel('Relative frequency');
% ylim([0 0.35])
% 
% 
% saveas(fig_ar,...
%     fullfile(save_path,[save_fname,'_ar.fig']),'fig');
% saveSameSize(fig_ar,'file',fullfile(save_path,[save_fname,'_ar.png']),...
%     'format','png','renderer','painters');
% epswrite(fullfile(save_path,[save_fname,'_ar.eps']));
