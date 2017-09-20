% 2015 11 24  Use real measurement az/el angles for model beampattern
% 2016 04 19  Use measured mean azimuth width to infer piston model diameter
%             Most processing and plotting codes are from rotate_all_click_20160419.m
% 2016 05 10  Re-process for plots in paper

clear
usrn = getenv('username');
if strcmp(usrn,'Wu-Jung')
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\EllipseDirectFit');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    addpath('F:\Dropbox\0_CODE\beampattern_processing');
else
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\EllipseDirectFit']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\beampattern_processing']);
end

% Set various paths
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

bat_proc_path = 'proc_output_rousettus_checked';
bat_proc_file_all = dir(fullfile(base_path,bat_proc_path,'rousettus_20150825_*'));
data_fname = 'all_model_rotated_3bats.mat';

azel_distr_path = 'az_el_distr_ellipse_20160419';
azel_distr_file = 'az_el_distr_ellipse_20160419_azel_xy_distr.mat';
azel_distr = load(fullfile(base_path,save_root,azel_distr_path,azel_distr_file));

A.base_path = base_path;
A.data_path = bat_proc_path;
A.save_path = save_path;
A.azel_distr_path = azel_distr_path;
A.azel_distr_file = azel_distr_file;

% Param
freq_model = 35e3;
it_shift_th = 0.005;
azel_bnd = [300,150];  % bound of which the rotate-shift difference larger than this would be discarded
cvec = 0:-3:-39;

plot_opt_all_click = 1;
plot_opt_indiv_click = 1;
save_opt = 1;

A.freq_modeled = freq_model;
A.iterative_shift_threshold = it_shift_th;
A.shift_bound_azel = azel_bnd;

% Set map projection
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);
[xxx,yyy] = mfwdtran(mstruct,[0,0,-90,90],[-180,180,0,0]);
xxx = xxx(1:2);
yyy = yyy(3:4);

A.map_projection = map_proj;
A.mstruct = mstruct;
A.map_plot_xlim = xxx;
A.map_plot_ylim = yyy;

% Set bp param
bp_info.c = 344;  % sound speed [m/s]
bp_info.freq = freq_model;  % [Hz]
bp_info.k = 2*pi*bp_info.freq/bp_info.c;  % wavenumber
bp_info.type = 'piston';

% Find aperture with best-fitting -3dB contour to mean azimuth of data
a_fine = (5:0.01:7)*1e-3; % [m]
theta_fine = 0:pi/1000:pi/2;
piston_bp = 20*log10(abs(2*besselj(1,bp_info.k*a_fine'*sin(theta_fine))./...
            (bp_info.k*a_fine'*sin(theta_fine))));
[C,~] = contour(theta_fine/pi*180,a_fine,piston_bp,[-3 -6],'fill','on');
c_level = parse_contour_output(C);
piston_3dB_c = [c_level(2).X; c_level(2).Y]';
[~,idx_mean] = min(abs(piston_3dB_c(:,1)-azel_distr.az_mean));
a_mean = piston_3dB_c(idx_mean,2);
close
bp_info.a = a_mean;

A.model_beampattern_info = bp_info;


% Model piston projection

% E_max_all = cell(length(bat_proc_file_all),1);
% ar_all = cell(length(bat_proc_file_all),1);
rotate_data = cell(length(bat_proc_file_all),1);
for iB = 1:length(bat_proc_file_all)

    bat_proc_file = bat_proc_file_all(iB).name;
    ss = strsplit(bat_proc_file,'_');
    save_fname = strjoin([script_name,ss(3:4)],'_');

    fprintf('File: %s\n',bat_proc_file);
    
    data = load(fullfile(base_path,bat_proc_path,bat_proc_file));
    good_call_idx = find(data.proc.chk_good_call);
    
    if plot_opt_all_click
        numrow = ceil(length(good_call_idx)/4);
        fig_clicks = figure;
        set(fig_clicks,'Position',[100 100 1000 140*numrow]);
    end
    
    for iC = good_call_idx'
 
        iC_save = find(iC==good_call_idx);

        % Get az/el from measurement
        [call_dB,az,el,ch_include_idx] = get_call_azel_dB_data(data,freq_model,iC);
        [~,mmidx] = max(call_dB);
        call_max_azel = [az(mmidx),el(mmidx)];
        
        % Model mic output
        piston_bp = model_beam(bp_info,call_max_azel,[az el]);
        
        % Fit ellipse
        [click_side,raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
            shift_rotate_bp(az,el,piston_bp','eckert4',it_shift_th,azel_bnd);
        
        % Check az, el shift and tilt
        if plot_opt_indiv_click==1
            fig_azel = figure;
            for iP=1:34
                azel = [raw.az(iP), raw.el(iP);...
                    rot_max.az(iP), rot_max.el(iP);...
                    rot_elpctr.az(iP), rot_elpctr.el(iP);...
                    rot_elpctr_tilt.az(iP), rot_elpctr_tilt.el(iP)];
                figure(fig_azel)
                plot(azel(:,1),azel(:,2),'.-');
                hold on
                text(azel(1,1),azel(1,2),num2str(iP));
                axis([-180 180 -90 90])
                grid on
                xlabel('Azimuth (deg)');
                ylabel('Elevation (deg)');
            end
            title(sprintf('%s, Call #%02d',regexprep(save_fname,'_','\\_'),iC));
            saveSameSize(fig_azel,'file',...
                fullfile(save_path,sprintf('%s_c%02d_azel_chk.png',save_fname,iC)),...
                'format','png','renderer','painters');
            close(fig_azel)
        end
        
        % Individual click shift/tilt contour & ellipse
        if plot_opt_indiv_click==1
            fig_fit = figure;
            set(fig_fit,'position',[140 100 900 480])
            
            for iSUB = 1:4
                switch iSUB
                    case 1
                        B = raw;
                        t_txt = 'Raw data';
                    case 2
                        B = rot_max;
                        t_txt = 'Best-fitting ellipse';
                    case 3
                        B = rot_elpctr;
                        t_txt = sprintf('Tilt %2.2fdeg',B.E.theta/pi*180);
                    case 4
                        B = rot_elpctr_tilt;
                        t_txt = 'Tilt compensated';
                end
                
                % valid vq_norm range after rotation
                if iSUB~=1
                    B.vq_norm(B.outofbnd_azel_idx_q) = NaN;
                end
                
                % bound for plotting ellipse
                if iSUB~=1
                    E = B.E;
                    c3db_xy = E.c3db_xy;
                    xmin = min(c3db_xy(:,1))-range(c3db_xy(:,1)*0.5);
                    xmax = max(c3db_xy(:,1))+range(c3db_xy(:,1)*0.5);
                    ymin = min(c3db_xy(:,2))-range(c3db_xy(:,1)*0.5);
                    ymax = max(c3db_xy(:,2))+range(c3db_xy(:,1)*0.5);
                end
                
                % plot contour and ellipse
                subplot(2,2,iSUB);
                contour(B.xq,B.yq,B.vq_norm,cvec(2:end),'fill','on');
                if iSUB~=1
                    hold on
                    fit_df = ezplot(E.eqt,[xmin,xmax,ymin,ymax]);  % plot ellipse
                    set(fit_df,'linecolor','b','linewidth',2);
                    plot(E.x0,E.y0,'r*');
                    text(E.x0,E.y0,sprintf('%2.3f, %2.3f',E.x0,E.y0))
                end
                axis equal; grid on
                xlim([xxx(1) xxx(2)])
                ylim([yyy(1) yyy(2)])
                title(t_txt);
%                 colorbar('ticks',fliplr(cvec(1:2:end)));
                colormap(parula(length(cvec)-1));
                caxis([cvec(end), 0])
            end
            
            suptitle(sprintf('%s, Call #%02d',regexprep(save_fname,'_','\\_'),iC));
            saveSameSize(fig_fit,'file',...
                fullfile(save_path,sprintf('%s_c%02d_contour.png',save_fname,iC)),...
                'format','png','renderer','painters');
            close(fig_fit)
        end

        % Rotated model output
        rotate_data{iB}(iC_save).raw = raw;
        rotate_data{iB}(iC_save).rot_max = rot_max;
        rotate_data{iB}(iC_save).rot_elpctr = rot_elpctr;
        rotate_data{iB}(iC_save).rot_elpctr_tilt = rot_elpctr_tilt;

        % Plot rotated bp using ellipse center
        if plot_opt_all_click
            figure(fig_clicks)
            subplot(numrow,4,find(iC==good_call_idx));
            contour(rot_elpctr_tilt.xq,rot_elpctr_tilt.yq,raw.vq_norm,cvec,'fill','on');
            axis equal
            grid on
            xlim([xxx(1) xxx(2)])
            ylim([yyy(1) yyy(2)])
            title(sprintf('Call #%02d',iC));
            colormap(parula(length(cvec)-1));
            caxis([cvec(end), 0])
        end

%         E_max_all{iF}(iC_save) = E_max;
%         ar_all{iF}(iC_save) = E_max.ar;
              
    end
    
    if plot_opt_all_click
        figure(fig_clicks)
        suptitle(sprintf('%s, all clicks',regexprep(save_fname,'_','\\_')));
        suptitle([save_fname,' all clicks']);
        saveSameSize(fig_clicks,'file',...
            fullfile(save_path,sprintf('%s_all_clicks.png',save_fname)),...
            'format','png','renderer','painters');
        close(fig_clicks)
    end

end

if save_opt==1
    A.rotate_data = rotate_data;
    save(fullfile(save_path,sprintf('%s_%s',script_name,data_fname)),'-struct','A');
end
