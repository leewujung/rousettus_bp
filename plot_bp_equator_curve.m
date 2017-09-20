% 2015 12 24  Plot bp curve along equator


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
    
    c3db_xy_all{iB} = cell(length(freq_wanted),length(good_call_idx));
    for iC = good_call_idx'

        iC_save = find(iC==good_call_idx);
        
        % Get call info and rotate beampattern using 35 kHz beampattern
        [call_dB,az,el,ch_include_idx] = get_call_azel_dB_data(data,35e3,iC);
        [click_side,raw,rot_max,rot_elp_ctr,rot_elp_ctr_tilt,fig_elp] = shift_rotate_bp(az(ch_include_idx),el(ch_include_idx),call_dB(ch_include_idx),'eckert4',0);
        
        % Get contour along equator
        [vq_xy,x,y] = get_equator_cut(rot_elp_ctr_tilt,call_dB);

