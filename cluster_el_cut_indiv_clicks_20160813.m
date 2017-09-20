% 2016 08 13  Cluster analysis of individual el-cut curves

clear

if isunix
    addpath('~/Dropbox/0_CODE/MATLAB/EllipseDirectFit');
    addpath('~/Dropbox/0_CODE/MATLAB/saveSameSize');
    
    data_base_path = '~/internal_2tb/Dropbox/0_ANALYSIS/bp_processing';
    save_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
else
    usrn = getenv('username');
    if strcmp(usrn,'Wu-Jung')
        addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
        addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    else
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
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

results_path = 'analysis_results_figs';
data_path = 'fig_el_cut_indiv_clicks_20160809';
knum = 2;
el_cut_all = [1,3,15,30];
freq_all = (25:10:55)*1e3;       % frequency to be plotted

for iF=1:length(freq_all)
    freq = freq_all(iF);
    
    for iE=1:length(el_cut_all)
        el_cut = el_cut_all(iE);
        
        
        data_file = sprintf('%s_elcut%02ddeg_%02dkHz.mat',...
            data_path,el_cut,freq/1e3);
        load(fullfile(data_base_path,results_path,data_path,data_file));
        
        [~,script_name,~] = fileparts(mfilename('fullpath'));
        save_path = fullfile(save_base_path,results_path,script_name);
        if ~exist(save_path,'dir')
            mkdir(save_path);
        end
        
        save_fname = sprintf('%s_knum%d_elcut%02d_freq%02dkHz',...
            script_name,knum,el_cut,freq/1e3);
        title_text = sprintf('elcut=%02ddeg, freq=%02dkHz, knum=%d',...
            el_cut,freq/1e3,knum);
        
        %  Clustering using curves between az=-60:60
        idx_left = kmeans(cut_az_left(:,41:81),knum);  
        idx_right = kmeans(cut_az_right(:,41:81),knum);
        
        cgrey = [1 1 1]*150/200;
        
        % Mean curves of all clusters
        fig_mean = figure('position',[50 150 1200 450]);
        for iK=1:knum
            subplot(121);
            plot(az_range,nanmean(cut_az_left(idx_left==iK,:),1),'linewidth',2);
            hold on
            subplot(122)
            plot(az_range,nanmean(cut_az_right(idx_right==iK,:),1),'linewidth',2);
            hold on
        end
        subplot(121); title('Left clicks');
        xlabel('Azimuth (deg)'); ylabel('Normalized beam energy');
        xlim([-180 180]); ylim([-30 3]); grid on;
        subplot(122); title('Right clicks');
        xlabel('Azimuth (deg)'); ylabel('Normalized beam energy');
        xlim([-180 180]); ylim([-30 3]); grid on;
        suptitle(title_text);
        
        % % saveas(fig_mean,fullfile(save_path,[save_fname,'_mean.fig']),'fig');
        saveSameSize(fig_mean,'file',fullfile(save_path,[save_fname,'_mean.png']),...
            'format','png','renderer','painters');
        % epswrite(fullfile(save_path,[save_fname,'_mean.eps']));
        
        % All curves in clusters
        fig_all = figure('units','normalized','outerposition',[0 0 1 1]);
        for iK=1:knum
            subplot(knum,2,1+2*(iK-1))
            plot(az_range,cut_az_left(idx_left==iK,:)','color',cgrey);
            hold on
            plot(az_range,nanmean(cut_az_left(idx_left==iK,:),1),'k','linewidth',2);
            xlim([-180 180]); ylim([-30 3]);
            if iK==1
                title('Left clicks');
            elseif iK==knum
                xlabel('Azimuth (deg)'); ylabel('Norm''ed energy');
            end
            
            subplot(knum,2,2+2*(iK-1))
            plot(az_range,cut_az_right(idx_right==iK,:)','color',cgrey);
            hold on
            plot(az_range,nanmean(cut_az_right(idx_right==iK,:),1),'k','linewidth',2);
            xlim([-180 180]); ylim([-30 3]);
            if iK==1
                title('Right clicks');
            elseif iK==knum
                xlabel('Azimuth (deg)'); ylabel('Norm''ed energy');
            end
        end
        suptitle(title_text);
        
%         saveas(fig_all,fullfile(save_path,[save_fname,'_all.fig']),'fig');
        saveSameSize_150(fig_all,'file',fullfile(save_path,[save_fname,'_all.png']),...
            'format','png','renderer','painters');
        % epswrite(fullfile(save_path,[save_fname,'_all.eps']));
        
        close(fig_all)
        close(fig_mean)
    end  % el_cut loop
end  % freq loop
