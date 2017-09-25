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
data_path = 'model_bp_proj_RaColony_multifreq_20170921_std1.0';
ss = strsplit(data_path,'_');
model_shape = ss{4};
noise_str = ss{end};


% Get data files
data_file = dir(fullfile(data_base_path,results_path,data_path,'*.mat'));

A.param.data_base_path = data_base_path;
A.param.data_path = data_path;
A.param.data_file = data_file;


% Load 1 file to get params
D = load(fullfile(data_base_path,results_path,data_path,data_file(1).name));  % rotated data
num_ch = size(D.v_mic,1);
[~,freqI] = sort(D.freq.all);  % order of sorted frequency
[xxx,yyy] = mfwdtran(D.map.mstruct,[0,0,-90,90],[-180,180,0,0]);
xy_lim = [xxx(1:2), yyy(3:4)];
num_freq = length(D.freq.all);
num_freq_plot = num_freq-2;

A.param.num_ch = num_ch;
A.map = D.map;


% Set plotting params
contour_sm_len = 10;
colorset = jet(num_freq_plot);
cgrey = 200*ones(1,3)/255;


% Set save folder
[~,script_name,~] = fileparts(mfilename('fullpath'));
script_name = sprintf('%s_%s',script_name,noise_str);
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

save_fname = script_name;



% Compile data from all clicks
for iS=1:length(data_file)
    fprintf('file %03d: %s\n',iS,data_file(iS).name);
    D = load(fullfile(data_base_path,results_path,...
                      data_path,data_file(iS).name));  % rotated data

    % Record click direction
    click_side(iS) = D.click_side;

    for iF=1:num_freq

        % NOTE: D.v_mic has already been normalized

        % interpolate for particular frequency --> rotated locations
        if isempty(D.rot_elpctr_tilt)
            continue
        end
        [~,vq_norm,azq,elq] = interp_bp(D.rot_elpctr_tilt.az/180*pi, ...
                                        D.rot_elpctr_tilt.el/180*pi,D.v_mic(:,iF),'rbf');
        azq = azq/pi*180;
        elq = elq/pi*180;

        % Get -3dB contour
        [azq_grid,elq_grid,vq_grid] = get_ortho_grid_azel(azq,elq,vq_norm);
        [~,c_level_nan] = get_main_contour(vq_grid,unique(azq_grid(:)), ...
                                                   unique(elq_grid(:)),-3);
        c3db_xy = [];
        [c3db_xy(:,1),c3db_xy(:,2)] = ...
            mfwdtran(D.map.mstruct,c_level_nan(:,2),c_level_nan(:,1));  % [az,el] to [x,y]
        c3db_xy_freq{iS,iF} = c3db_xy;
        
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
            shift_rotate_bp(D.rot_elpctr_tilt.az/180*pi,...
                            D.rot_elpctr_tilt.el/180*pi,...
                            D.v_mic(:,iF),D.map.map_projection,0.005);
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
end  % all click files

A.bpctr = bpctr;
A.c3db_xy_freq = c3db_xy_freq;
A.click_side = click_side;
A.rot_n = rot_n;


% Save results
save_fname = [script_name,'_results.mat'];
save(fullfile(save_path,save_fname),'-struct','A');



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
    fig_bpctr = figure('position',[100,100,1100,700]);

    for ii=1:4
        subplot(2,2,ii)
        axesm(D.map.map_projection);
        gridm('gcolor',cgrey,'glinestyle','-');
        framem('fedgecolor',cgrey,'flonlimit',[-180 180]);
        axis off
        hold on
    end

    cnt = 0;
    for iF=freqI(2:num_freq-1)  % 25:5:55 kHz, need to use sorted freq sequence
        cnt = cnt+1;
        for iS=1:length(data_file)
            % bpctr
            if click_side(iS)==1
                subplot(222);
            else
                subplot(221);
            end
            if ~isempty(rot_n(iS,iF))
                switch bpctr_opt
                  case 'max'
                    plotm(bpctr(iS,iF).max_el,bpctr(iS,iF).max_az,'.','markersize',8,...
                          'linewidth',2,'color',colorset(cnt,:));
                  case 'top'
                    plotm(bpctr(iS,iF).top_el,bpctr(iS,iF).top_az,'.','markersize',8,...
                          'linewidth',2,'color',colorset(cnt,:));
                  case 'ectr'
                    plotm(bpctr(iS,iF).ectr_el,bpctr(iS,iF).ectr_az,'.','markersize',8,...
                          'linewidth',2,'color',colorset(cnt,:));
                end
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
            plot(xy_sm(:,1),xy_sm(:,2),'linewidth',0.5,'color',colorset(cnt,:));
            clear xy_sm
        end
    end
    for ii=1:4
        subplot(2,2,ii);
        tightmap
        colormap(jet(num_freq_plot))
        colorbar('Ticks',linspace(0+1/num_freq_plot/2,1-1/num_freq_plot/2,num_freq_plot),...
                 'TickLabels',{num2str(D.freq.all(freqI(2:end-1))'/1e3)},'location','southoutside');
        grid
    end


    % Save figures
    saveas(fig_bpctr,fullfile(save_path,sprintf('%s_%s.fig',save_fname,bpctr_opt)),'fig');
    saveSameSize_res(fig_bpctr,150,'file',fullfile(save_path,sprintf('%s_%s.png',save_fname,bpctr_opt)),...
                     'format','png','renderer','painters');
    epswrite(fullfile(save_path,sprintf('%s_%s.eps',save_fname,bpctr_opt)));

end