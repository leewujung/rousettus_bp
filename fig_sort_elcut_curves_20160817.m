% 2015 12 27  Make multi-freq composite clicks
% 2016 04 20  Update to make it compatible with new rotated click format
% 2016 05 07  Reprocess for plots for NIFTI poster
% 2016 08 08  Update rotated data format after 20160721
% 2016 08 13  Composite click according to clustering results

clear

usrn = getenv('username');
if strcmp(usrn,'Wu-Jung')
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
else
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
end

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    bat_proc_path = 'F:\Dropbox\0_ANALYSIS\bp_processing\proc_output_rousettus_checked';
    data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
else
    bat_proc_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\proc_output_rousettus_checked'];
    data_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
    save_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
end
results_path = 'analysis_results_figs';

% Set save folder
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Load el-cut curves
bat_num = 'all';  % all,34271,36134,39184
el_cut = 3;
freq_cluster = 35e3;
elcut_path = 'fig_el_cut_indiv_clicks_20160818';
elcut_file = sprintf('%s_bat%s_elcut%02ddeg_%02dkHz.mat',...
    elcut_path,bat_num,el_cut,freq_cluster/1e3);
C = load(fullfile(data_base_path,results_path,elcut_path,elcut_file));

save_fname = sprintf('%s_bat%s_elcut%02ddeg_sort',...
    script_name,bat_num,el_cut);

% Sort the curves according to slope
th = -6;
vert_shift = 2;
m_fac = 1;
cvec = 0:-1:-30;
fig_img = figure('position',[220 280 1050 430]);
fig_cntr = figure('position',[220 280 1050 430]);
fig_wtrf = figure('position',[220 280 1050 430]);
fig_lines = figure('position',[610 120 780 820]);
for iSide=1:2
    if iSide==1 % left
%         sort_idx_range = 41:65;
        sort_idx_range = 47:65;
        idx_want = ~isnan(C.cut_az_left(:,sort_idx_range(1)));
        cut_el_want = C.cut_az_left(idx_want,:);
        th = th;
    else  % right
        sort_idx_range = 57:75;
        idx_want = ~isnan(C.cut_az_right(:,sort_idx_range(end)));
        cut_el_want = C.cut_az_right(idx_want,:);
        th = th;
    end
    [~,min_idx] = min(abs(cut_el_want(:,sort_idx_range)-th),[],2);
    min_idx = min_idx+sort_idx_range(1)-1;
    [min_sort,min_sort_idx] = sort(min_idx);
    num_click = length(min_sort_idx);
    num_click_round = ceil(num_click/10)*10;
    
    figure(fig_img)
    subplot(1,2,iSide);
    if iSide==2
        imagesc(C.az_range,1:length(min_sort_idx),...
            flipud(cut_el_want(min_sort_idx,:)));
        hold on
        plot(C.az_range(min_idx(min_sort_idx)),...
            fliplr(1:length(min_sort_idx)),...
            'w','linewidth',2);
    else
        imagesc(C.az_range,1:length(min_sort_idx),...
            cut_el_want(min_sort_idx,:));
        hold on
        plot(C.az_range(min_idx(min_sort_idx)),...
            1:length(min_sort_idx),...
            'w','linewidth',2);
    end
    caxis([-40 0])
    colorbar

    figure(fig_cntr)
    subplot(1,2,iSide);
    if iSide==2
        contourf(C.az_range,1:num_click,...
            cut_el_want(min_sort_idx,:),cvec);
        hold on
        plot(C.az_range(min_idx(min_sort_idx)),...
            1:num_click,...
            'w','linewidth',2);
    else
        contourf(C.az_range,1:num_click,...
            cut_el_want(flipud(min_sort_idx),:),cvec);
        hold on
        plot(C.az_range(min_idx(min_sort_idx)),...
            fliplr(1:num_click),...
            'w','linewidth',2);
    end
    colormap(parula(length(cvec)));
    colorbar
    
    figure(fig_wtrf)
    subplot(1,2,iSide);
    if iSide==2
        waterfall(C.az_range,1:num_click,...
            cut_el_want(min_sort_idx,:));
        view([20 55])
    else
        waterfall(C.az_range,1:num_click,...
            cut_el_want(flipud(min_sort_idx),:));
        view([-20 55])
    end
    xlabel('Azimuth (deg)')
    ylabel('Click number')
    zlabel('Normalized beam energy')
    caxis([-40 0])
    zlim([-25 5])
    ylim([0 num_click_round])
    xlim([-180 180])
    colorbar

    figure(fig_lines)
    subplot(1,2,iSide);
    if iSide==2
        plot(C.az_range,...
            cut_el_want(min_sort_idx,:)*m_fac+...
            repmat(vert_shift*(0:num_click-1)',1,size(cut_el_want,2)),...
            'k');
        ylim([-20 310])
    else
        plot(C.az_range,...
            cut_el_want(flipud(min_sort_idx),:)*m_fac+...
            repmat(vert_shift*(0:num_click-1)',1,size(cut_el_want,2)),...
            'k');
        ylim([-20 260])
    end
    xlim([-180 180])
%     ylim([-20 vert_shift*num_click_round])
    xlabel('Azimuth (deg)');
    ylabel('Click number (sorted)');
    set(gca,'ytick',(0:vert_shift*20:vert_shift*(num_click_round+10))-20,...
        'yticklabel',0:20:num_click_round+10,...
        'xtick',-180:60:180);
end

title_text = sprintf('bat %s, elcut %d deg, th %d dB',...
    bat_num,el_cut,th);

figure(fig_img);
mtit(title_text);
saveSameSize(fig_img,'file',fullfile(save_path,[save_fname,'_imgsc.png']),...
    'format','png','renderer','painters');
% epswrite(fullfile(save_path,[save_fname,'_imgsc.eps']));

figure(fig_cntr);
mtit(title_text);
saveSameSize(fig_cntr,'file',fullfile(save_path,[save_fname,'_cntr.png']),...
    'format','png','renderer','painters');
% epswrite(fullfile(save_path,[save_fname,'_cntr.eps']));

figure(fig_wtrf);
mtit(title_text);
saveSameSize(fig_wtrf,'file',fullfile(save_path,[save_fname,'_wtrf.png']),...
    'format','png','renderer','painters');
epswrite(fullfile(save_path,[save_fname,'_wtrf.eps']));

figure(fig_lines);
mtit(title_text);
saveSameSize(fig_lines,'file',fullfile(save_path,[save_fname,'_lines.png']),...
    'format','png','renderer','painters');
epswrite(fullfile(save_path,[save_fname,'_cntr.eps']));
