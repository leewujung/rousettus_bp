% 2015 12 28  Plots for multi-freq composite clicks
% 2016 04 20  Separate this part out from multifreq_composite_click.m
% 2016 05 07  Plot for NIFTI poster
% 2016 05 07  Plot for ASA talk
% 2016 07 23  Plot for paper; best-fitting ellipse for composite click

clear

if isunix
else
    usrn = getenv('username');
    if strcmp(usrn,'Wu-Jung')
        addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    else
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
    end
    
    if strcmp(usrn,'Wu-Jung')
        data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
        save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    else
        data_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
        save_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
    end
end

save_plot_opt = 1;
results_path = 'analysis_results_figs';
data_path = 'multifreq_composite_click_20160420';
it_shift_th = 0.005;
azel_bnd = [300,150];

% Save path
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(data_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

S.param.data_path = data_path;
S.param.it_shift_th = it_shift_th;
S.param.azel_bnd = azel_bnd;

% Each individual bat
bat = {'3bat','36134','34271','39184'};
for iBAT=1%:length(bat)
    data_file = ['multifreq_composite_click_20160420_all_clicks_',bat{iBAT},'.mat'];
    D = load(fullfile(data_base_path,results_path,data_path,data_file));
    save_data_file = sprintf('%s_results.mat',script_name);
    
    S.param.data_file = data_file;
    S.param.freq = D.param.freq_wanted;
    
    num_freq = length(D.param.freq_wanted);
    for iF=1:num_freq
        freq_wanted = D.param.freq_wanted(iF);
        
        % params for figure
        vq_norm_min = -27;
        cvec = 0:-3:(floor(vq_norm_min/3)-1)*3;
        cvec_min_idx = find(cvec-vq_norm_min<0,1,'first');
        map_proj = D.param.map.map_projection;
        
        S.cvec = cvec;
        S.param.map = D.param.map;

        % Right click =============
        title_text = sprintf('Composite click, %d kHz, right click',freq_wanted/1e3);
        save_fname = sprintf('%s_%dkHz_right',script_name,freq_wanted/1e3);

        azq = D.averaged_composite.right.interp(iF).azq_avg;
        elq = D.averaged_composite.right.interp(iF).elq_avg;
        vq_norm = D.averaged_composite.right.interp(iF).vq_norm_avg;
        
        [raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
            shift_rotate_bp_composite(azq,elq,vq_norm,map_proj,it_shift_th,azel_bnd);
        S.right(iF).raw = raw;
        S.right(iF).rot_max = rot_max;
        S.right(iF).rot_elpctr = rot_elpctr;
        S.right(iF).rot_elpctr_tilt = rot_elpctr_tilt;

        % Plot rotation procedure
        fig_fit = plot_indiv_click_rotate(S.right(iF),cvec,D.param.map.mstruct);
        suptitle(title_text);
        saveSameSize(fig_fit,'file',...
            fullfile(save_path,save_fname),...
            'format','png','renderer','painters');
        close(fig_fit)

        % Left click =============
        title_text = sprintf('Composite click, %d kHz, left click',freq_wanted/1e3);
        save_fname = sprintf('%s_%dkHz_left',script_name,freq_wanted/1e3);

        azq = D.averaged_composite.left.interp(iF).azq_avg;
        elq = D.averaged_composite.left.interp(iF).elq_avg;
        vq_norm = D.averaged_composite.left.interp(iF).vq_norm_avg;
        
        [raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
            shift_rotate_bp_composite(azq,elq,vq_norm,map_proj,it_shift_th,azel_bnd);
        S.left(iF).raw = raw;
        S.left(iF).rot_max = rot_max;
        S.left(iF).rot_elpctr = rot_elpctr;
        S.left(iF).rot_elpctr_tilt = rot_elpctr_tilt;

        % Plot rotation procedure
        fig_fit = plot_indiv_click_rotate(S.left(iF),cvec,D.param.map.mstruct);
        suptitle(title_text);
        saveSameSize(fig_fit,'file',...
            fullfile(save_path,save_fname),...
            'format','png','renderer','painters');
        close(fig_fit)
    end
end

save(fullfile(save_path,save_data_file),'-struct','S');


% Plot best-fitting ellipse param
for iF=1:num_freq
    [ar_left(iF),xy] = get_ellipse_ar(S.left(iF).rot_elpctr_tilt.E);
    x_max = range(xy(:,1))/2;
    y_max = range(xy(:,2))/2;
    [el,az] = minvtran(S.param.map.mstruct,x_max,y_max);
    azel_ar_left(iF) = el/az;

    [ar_right(iF),xy] = get_ellipse_ar(S.right(iF).rot_elpctr_tilt.E);
    x_max = range(xy(:,1))/2;
    y_max = range(xy(:,2))/2;
    [el,az] = minvtran(S.param.map.mstruct,x_max,y_max);
    azel_ar_right(iF) = el/az;
end

figure
plot(S.param.freq/1e3,ar_left,'-o');
hold on
plot(S.param.freq/1e3,ar_right,'-o');
legend('Left','Right')
grid
xlabel('Frequency (kHz)');
ylabel('XY aspect ratio');
title('XY aspect ratio across frequency');
saveSameSize(gcf,'file',...
    fullfile(save_path,sprintf('%s_xy_ar.png',script_name)),...
    'format','png','renderer','painters');

figure
plot(S.param.freq/1e3,azel_ar_left,'-o');
hold on
plot(S.param.freq/1e3,azel_ar_right,'-o');
legend('Left','Right')
grid
xlabel('Frequency (kHz)');
ylabel('AZ/EL aspect ratio');
title('AZ/EL aspect ratio across frequency');
saveSameSize(gcf,'file',...
    fullfile(save_path,sprintf('%s_azel_ar.png',script_name)),...
    'format','png','renderer','painters');





