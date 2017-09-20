% 2015 12 03  Merge data from all clicks
% 2015 12 12  Use 35 kHz rotation for other frequencies
% 2015 12 22  Adopt the merging click code for finding multi-freq center of
%             the beampatterns
% 2015 12 24  Move the plotting part to 'compile_info_plot.m'
%             This is only for processing

clear
warning off

usrn = getenv('username');
base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
save_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20151222_multifreq_bp_center'];
if ~exist(save_path,'dir')
    mkdir(save_path);
end
bat_proc_path = './proc_output_rousettus_checked';
bat_proc_file_all = dir(fullfile(base_path,bat_proc_path,'rousettus_20150825_36134*'));
plot_opt = 0;
freq_wanted = [35,20:5:30,40:5:50]*1e3;
num_freq = length(freq_wanted);

addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\EllipseDirectFit']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);


% Set map projection
mstruct = defaultm('eckert4');
mstruct = defaultm(mstruct);
        
% Process all files
az_ectr_all = cell(length(bat_proc_file_all),num_freq);
el_ectr_all = az_ectr_all;
az_max_all = az_ectr_all;
el_max_all = az_ectr_all;
az_ectr_tilt_all = az_ectr_all;
el_ectr_tilt_all = az_ectr_all;
x_ectr_all = az_ectr_all;
y_ectr_all = az_ectr_all;
call_dB_norm_all = az_ectr_all;
click_side_all = az_ectr_all;  % 1/0 for each channel in each click
click_side_single_all = cell(length(bat_proc_file_all),1);  % 1/0 for each click
elp_ctr_x_all = cell(length(bat_proc_file_all),1);
elp_ctr_y_all = cell(length(bat_proc_file_all),1);
c3db_xy_all = cell(length(bat_proc_file_all),1);
equator_cut_curve_all = cell(length(bat_proc_file_all),1);
equator_cut_x_all = cell(length(bat_proc_file_all),1);
equator_cut_y_all = cell(length(bat_proc_file_all),1);
for iB = 1:length(bat_proc_file_all)

    bat_proc_file = bat_proc_file_all(iB).name;
    fprintf('File: %s\n',bat_proc_file);
    
    data = load(fullfile(base_path,bat_proc_path,bat_proc_file));
    mic_num = 1:data.mic_data.num_ch_in_file;
    good_call_idx = find(data.proc.chk_good_call);
    
    if plot_opt
        fig_all = figure;
    end
    numrow = ceil(length(good_call_idx)/2);
    
    equator_cut_curve_all{iB} = cell(length(good_call_idx),1);
    equator_cut_x_all{iB} = cell(length(good_call_idx),1);
    equator_cut_y_all{iB} = cell(length(good_call_idx),1);
    c3db_xy_all{iB} = cell(length(freq_wanted),length(good_call_idx));
    for iC = good_call_idx'

        iC_save = find(iC==good_call_idx);
        
        % Get call info and rotate beampattern using 35 kHz beampattern
        [call_dB,az,el,ch_include_idx] = get_call_azel_dB_data(data,35e3,iC);
        [click_side,raw,rot_max,rot_elp_ctr,rot_elp_ctr_tilt,fig_elp] = shift_rotate_bp(az(ch_include_idx),el(ch_include_idx),call_dB(ch_include_idx),'eckert4',0);
        
        % Get contour along equator
        [vq_xy,x,y] = get_equator_cut(rot_elp_ctr_tilt,call_dB,mstruct);
        equator_cut_curve_all{iB}{iC_save} = vq_xy;
        equator_cut_x_all{iB}{iC_save} = x;
        equator_cut_y_all{iB}{iC_save} = y;

%         figure;
%         surf(rot_elp_ctr_tilt.xq,rot_elp_ctr_tilt.yq,raw.vq_norm,'edgecolor','none');
%         hold on
%         plot3(x,y,vq_xy,'.');
         
        % Save data
        az_max_all{iB}(iC_save,:) = rot_max.az;
        el_max_all{iB}(iC_save,:) = rot_max.el;
        
        az_ectr_all{iB}(iC_save,:) = rot_elp_ctr.az;
        el_ectr_all{iB}(iC_save,:) = rot_elp_ctr.el;
        x_ectr_all{iB}(iC_save,:) = rot_elp_ctr.x;
        y_ectr_all{iB}(iC_save,:) = rot_elp_ctr.y;
        
        az_ectr_tilt_all{iB}(iC_save,:) = rot_elp_ctr_tilt.az;
        el_ectr_tilt_all{iB}(iC_save,:) = rot_elp_ctr_tilt.el;

        % Multi-freq analysis
        az_ecen_tilt = az_ectr_tilt_all{iB}(iC_save,:);
        el_ecen_tilt = el_ectr_tilt_all{iB}(iC_save,:);

        elp_ctr_x_local = nan(num_freq,1);
        elp_ctr_y_local = nan(num_freq,1);
        c3db_avg_x = nan(num_freq,1);
        c3db_avg_y = nan(num_freq,1);

        fig_elp_f = figure;
        for iF=1:num_freq
            [call_dB_f,az_f,el_f,ch_include_idx_f] = get_call_azel_dB_data(data,freq_wanted(iF),iC);  % get good mic index
            call_dB_norm_all{iB,iF}(iC_save,:) = call_dB_f - max(call_dB_f);  % save call_dB_norm data
            ch_include_idx_f = ch_include_idx_f & ~isnan(az_ecen_tilt);
            [~,vq_norm_f,azq_f,elq_f] = interp_bp(az_ecen_tilt(ch_include_idx_f)/180*pi,el_ecen_tilt(ch_include_idx_f)/180*pi,call_dB_f(ch_include_idx_f),'rbf');  % use the first frequency data for finding ellipse center
            [xq_f,yq_f] = mfwdtran(mstruct,elq_f/pi*180,azq_f/pi*180);
            
            % center of the best-fitting ellipse
            figure(fig_elp_f);
            E_max = plot_bp_fit_ellipse(gca,xq_f,yq_f,vq_norm_f);
            elp_ctr_x_local(iF,:) = E_max.x0;
            elp_ctr_y_local(iF,:) = E_max.y0;
            
            % Exrtract -3dB contour
            [~,c_main_nan] = get_main_contour(vq_norm_f,azq_f(1,:)/pi*180,elq_f(:,1)/pi*180,-3);
            [c3db_xy(:,1),c3db_xy(:,2)] = mfwdtran(mstruct,c_main_nan(:,2),c_main_nan(:,1));
            c3db_xy_all{iB}{iF,iC_save} = c3db_xy;

            % Mean point of -3db contour
            c3db_avg_x(iF,:) = mean(c3db_xy(:,1));
            c3db_avg_y(iF,:) = mean(c3db_xy(:,2));
            
            clear c3db_xy
        end
        close(fig_elp_f);
        
        click_side_single_all{iB}(iC_save) = click_side;
        click_side_all{iB}(iC_save,:) = click_side*ones(length(az),1);
        elp_ctr_x_all{iB}(:,iC_save) = elp_ctr_x_local;
        elp_ctr_y_all{iB}(:,iC_save) = elp_ctr_y_local;
        
        % Plot rotated bp using ellipse center
        if plot_opt
            figure(fig_all)
            subplot(2,numrow,find(iC==good_call_idx));
            [C,~] = contour(rot_elp_ctr_tilt.xq,rot_elp_ctr_tilt.yq,raw.vq_norm,0:-3:-39,'fill','on');
            title(sprintf('Call #%02d',iC));
        end
    end
    
    if plot_opt
        figure(fig_all)
        suptitle(sprintf('%s call #%02d',bat_proc_file,iC));
    end
    
end
warning on

save(fullfile(save_path,'all_click_data.mat'));

