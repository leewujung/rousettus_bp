% 2017 09 24  Use rotate/shift output to gather stats of beam center for all frequencies

clear

usrn = getenv('username');
if isunix
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/rbfinterp_v1.2');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/saveSameSize');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/EllipseDirectFit');
    bat_proc_path = '~/internal_2tb/Dropbox/0_ANALYSIS/bp_processing/proc_output_rousettus_new_checked';
    data_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
    save_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
else
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    bat_proc_path = 'F:\Dropbox\0_ANALYSIS\bp_processing\proc_output_rousettus_new_checked';
    data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
end

results_path = 'analysis_results_figs';
data_path = 'rotate_all_click_2010308';
bat_proc_fpre = 'rousettus_20150825';


% Set params
freq_wanted = (20:5:60)*1e3;
bat_num = 'all';  % all,34271,36134,39184
num_freq = length(freq_wanted);
num_freq_plot = num_freq-2;
contour_sm_len = 10;
colorset = jet(num_freq-2);
cgrey = 200*ones(1,3)/255;

A.param.freq_wanted = freq_wanted;
A.param.bat_num = bat_num;


% Get data files
if strcmp(bat_num,'all')
    data_file = dir(fullfile(data_base_path,results_path,data_path,'*.mat'));
else
    data_file = dir(fullfile(data_base_path,results_path,data_path,['*',bat_num,'*.mat']));
end

A.param.data_base_path = data_base_path;
A.param.rotated_data_path = data_path;
A.param.bat_proc_path = bat_proc_path;
A.param.rotated_data_file = data_file;


% Load 1 file to get params
D = load(fullfile(data_base_path,results_path,data_path,data_file(1).name));  % rotated data
num_ch = length(D.raw_meas.call_dB);
[xxx,yyy] = mfwdtran(D.map.mstruct,[0,0,-90,90],[-180,180,0,0]);
xy_lim = [xxx(1:2), yyy(3:4)];

A.param.num_ch = num_ch;
A.map = D.map;


% Set save folder
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

save_fname = sprintf('%s_bat%s',script_name,bat_num);


% Compile data from all clicks
for iS=1:length(data_file)
    fprintf('file %03d: %s\n',iS,data_file(iS).name);
    D = load(fullfile(data_base_path,results_path,...
                      data_path,data_file(iS).name));  % rotated data

    % Record click direction
    if isempty(D.raw_meas.click_side) || isnan(D.raw_meas.click_side)
        click_side(iS) = NaN;
    else
        click_side(iS) = D.raw_meas.click_side;
    end

    % Load multi-freq data
    ss = strsplit(strtok(data_file(iS).name,'.'),'_');
    bat_num = ss{5};
    trial_num = ss{6};
    click_num = str2double(ss{7}(2:3));
    bat_proc_fname = sprintf('%s_%s_%s_mic_data_bp_proc.mat',...
                             bat_proc_fpre,bat_num,trial_num);
    data = load(fullfile(bat_proc_path,bat_proc_fname));
    fprintf('--- corresponding original data file: %s\n',bat_proc_fname);
    for iF=1:num_freq

        [call_dB,az,el,ch_include_idx] = ...
            get_call_azel_dB_data(data,freq_wanted(iF),click_num);  % get good mic index
        call_dB = call_dB - max(call_dB);
        call_dB(~ch_include_idx) = NaN;
        call_dB_norm = call_dB - max(call_dB);  % normalize

        % interpolate for particular frequency --> rotated locations
        if isempty(D.rot_elpctr_tilt)
            c3db_xy_freq{iS,iF} = [NaN,NaN];
            bpctr(iS,iF).max_el = NaN;
            bpctr(iS,iF).max_az = NaN;
            bpctr(iS,iF).top_el = NaN;
            bpctr(iS,iF).top_az = NaN;
            bpctr(iS,iF).ectr_el = NaN;
            bpctr(iS,iF).ectr_az = NaN;
            rot_n(iS,iF) = NaN;
            continue;
        end
        [~,vq_norm,azq,elq] = interp_bp(D.rot_elpctr_tilt.az/180*pi, ...
                                        D.rot_elpctr_tilt.el/180*pi,call_dB,'rbf');
        azq = azq/pi*180;
        elq = elq/pi*180;

        % Get -3dB contour
        [azq_grid,elq_grid,vq_grid] = get_ortho_grid_azel(azq,elq,vq_norm);
        [~,c_level_nan] = get_main_contour(vq_grid,unique(azq_grid(:)),unique(elq_grid(:)),-3);
        [c3db_xy(:,1),c3db_xy(:,2)] = ...
            mfwdtran(D.map.mstruct,c_level_nan(:,2),c_level_nan(:,1));  % [az,el] to [x,y]

        c3db_xy_freq{iS,iF} = c3db_xy;
        clear c3db_xy
        
        % Find beam center --> rotated locations (using azq/elq from above)
        % --- max beam energy point
        xx = vq_norm(:);
        [~,vq_norm_max_idx] = max(xx);
        bpctr(iS,iF).max_el = elq(vq_norm_max_idx);
        bpctr(iS,iF).max_az = azq(vq_norm_max_idx);
        % --- averaged loc of all >-1 dB points
        xx(isnan(xx)) = -inf;
        [~,sort_idx] = sort(xx,'descend');
        ii = xx(sort_idx)>-1;
        bpctr(iS,iF).top_el = mean(elq(sort_idx(ii)));
        bpctr(iS,iF).top_az = mean(azq(sort_idx(ii)));
        % --- center of best-fitting ellipse
        [~,raw,rot_max,rot_elpctr,rot_elpctr_tilt,rot_n(iS,iF)] = ...
            shift_rotate_bp(D.rot_elpctr_tilt.az(ch_include_idx)/180*pi,...
                            D.rot_elpctr_tilt.el(ch_include_idx)/180*pi,...
                            call_dB(ch_include_idx),D.map.map_projection,0.005);
        if isempty(rot_max)
            bpctr(iS,iF).ectr_el = NaN;
            bpctr(iS,iF).ectr_az = NaN;
        else
            if rot_max.E.x0<xy_lim(1) || rot_max.E.x0>xy_lim(2) ||...
                    rot_max.E.y0<xy_lim(3) || rot_max.E.y0>xy_lim(4)
                bpctr(iS,iF).ectr_el = [];
                bpctr(iS,iF).ectr_az = [];
            else            
                [el_ectr,az_ectr] = minvtran(D.map.mstruct,rot_max.E.x0,rot_max.E.y0);  % inverse map projection
                [bpctr(iS,iF).ectr_el,bpctr(iS,iF).ectr_az] = rotatem(el_ectr,az_ectr,...
                                     [bpctr(iS,iF).max_el,bpctr(iS,iF).max_az],...
                                                                  'inverse','degrees');
            end
        end
    end  % all freq

    if mod(iS,10)==0  % save results everything 10 files
        A.bpctr = bpctr;
        A.c3db_xy_freq = c3db_xy_freq;
        A.click_side = click_side;
        A.rot_n = rot_n;
        save(fullfile(save_path,sprintf('%s_file001-%03d_.mat',save_fname,iS)),...
             '-struct','A');
    end
end  % all click files

A.bpctr = bpctr;
A.c3db_xy_freq = c3db_xy_freq;
A.click_side = click_side;
A.rot_n = rot_n;



% Save results
save(fullfile(save_path,[save_fname,'_results.mat']),'-struct','A');



% Plot all bpctr and contours
for iB=1:3
    switch iB
      case 1
        bpctr_opt = 'max';
      case 2
        bpctr_opt = 'top';
      case 3
        bpctr_opt = 'ectr';
    end
    num_freq_plot = num_freq-2;
    fig_bpctr = figure('position',[100,100,1100,700]);

    for ii=1:4
        subplot(2,2,ii)
        axesm(D.map.map_projection);
        gridm('gcolor',cgrey,'glinestyle','-');
        framem('fedgecolor',cgrey,'flonlimit',[-180 180]);
        axis off
        hold on
    end

    for iF=2:num_freq-1  % 25:5:55 kHz
        for iS=1:length(data_file)
            % bpctr
            if click_side(iS)==1
                subplot(222);
            else
                subplot(221);
            end
            switch bpctr_opt
              case 'max'
                plotm(bpctr(iS,iF).max_el,bpctr(iS,iF).max_az,'.','markersize',8,...
                      'linewidth',2,'color',colorset(iF-1,:));
              case 'top'
                plotm(bpctr(iS,iF).top_el,bpctr(iS,iF).top_az,'.','markersize',8,...
                      'linewidth',2,'color',colorset(iF-1,:));
              case 'ectr'
                plotm(bpctr(iS,iF).ectr_el,bpctr(iS,iF).ectr_az,'.','markersize',8,...
                      'linewidth',2,'color',colorset(iF-1,:));
            end
            % -3dB contours
            if click_side(iS)==1
                subplot(224);
            else
                subplot(223);
            end
            xy = c3db_xy_freq{iS,iF};
            xy_sm(:,1) = smooth(xy(:,1),contour_sm_len);
            xy_sm(:,2) = smooth(xy(:,2),contour_sm_len);
            xy_sm(isnan(xy(:,1)),:) = NaN;
            plot(xy_sm(:,1),xy_sm(:,2),'linewidth',0.5,'color',colorset(iF-1,:));
            clear xy_sm
        end
    end

    for ii=1:4
        subplot(2,2,ii);
        tightmap
        colormap(jet(num_freq_plot))
        colorbar('Ticks',linspace(0+1/num_freq_plot/2,1-1/num_freq_plot/2,num_freq_plot),...
                 'TickLabels',{num2str(freq_wanted(2:end-1)'/1e3)},'location','southoutside');
        grid
    end


    % Save figures
    saveas(fig_bpctr,fullfile(save_path,sprintf('%s_%s.fig',save_fname,bpctr_opt)),'fig');
    saveSameSize_res(fig_bpctr,150,'file',fullfile(save_path,sprintf('%s_%s.png',save_fname,bpctr_opt)),...
                     'format','png','renderer','painters');
    epswrite(fullfile(save_path,sprintf('%s_%s.eps',save_fname,bpctr_opt)));

end