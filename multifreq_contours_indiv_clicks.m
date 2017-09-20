% 2015 12 27  Plot contours of the same frequency across all clicks together

clear

usrn = getenv('username');
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);

plot_opt = 1;
save_opt = 1;

% Load compiled rotated data
base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
compile_path = '20151224_all_click_rotate';
compile_file = 'all_click_rotated_39184.mat';
load(fullfile(base_path,compile_path,compile_file));

% Fix path with current username
if ~strcmp(usrn,'wlee76')
    if strfind(base_path, 'wlee76')
        base_path = regexprep(base_path,'wlee76',usrn);
    end
end

A.base_path = base_path;
A.compile_path = compile_path;
A.compile_file = compile_file;

% Path/filename for saving data and figures
save_header = 'rousettus_39184';
save_path = fullfile(base_path,'20151228_indiv_click_contour');  % set path for saving files
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Frequency to be plotted
freq_wanted = (20:5:50)*1e3;
num_freq = length(freq_wanted);

A.freq_wanted  = freq_wanted;


warning off
% Get contours from all freq
c3db_xy_all = cell(length(processed_files),1);
for iB=1:length(processed_files)
    bat_proc_file = processed_files(iB).name;
    ss = strsplit(bat_proc_file,'_');
    save_fname = strjoin(ss(3:4),'-');
    fprintf('File: %s\n',bat_proc_file);
    
    data = load(fullfile(base_path,data_path,bat_proc_file));
    good_call_idx = find(data.proc.chk_good_call);

    for iC = good_call_idx'
        iC_save = find(iC==good_call_idx);
        fprintf('Call %02d\n',iC);
        
        for iF=1:num_freq
            % Get data and good mic index
            [call_dB,az,el,ch_include_idx] = get_call_azel_dB_data(data,freq_wanted(iF),iC);
            call_dB_norm = call_dB - max(call_dB);  % normalize
            [~,vq_norm,~,~] = interp_bp(az,el,call_dB,'rbf');  % use the first frequency data for finding ellipse center
            
            % Get -3dB contour
            M = rotate_data{iB}(iC_save).rot_elpctr_tilt;
            idx_notnan = ~isnan(M.azq);  % index of non-NaN data
            azq = M.azq(idx_notnan);
            elq = M.elq(idx_notnan);
            azq = min(azq(:)):max(azq(:));
            elq = min(elq(:)):max(elq(:));
            [azq,elq] = meshgrid(azq,elq);
            [azq,elq,vq] = griddata(M.azq(idx_notnan),M.elq(idx_notnan),vq_norm(idx_notnan),azq,elq);  % griddata to create regular orthogonal azq/elq
            
            [~,c_main_nan] = get_main_contour(vq,unique(azq),unique(elq),-3);  % get main contour at -3dB with NaN insert for break contour
            [c3db_xy(:,1),c3db_xy(:,2)] = mfwdtran(mstruct,c_main_nan(:,2),c_main_nan(:,1));  % [az,el] to [x,y]
            c3db_xy_all{iB}{iF,iC_save} = c3db_xy;
                        
            clear c3db_xy
        end
    end
end
A.contour_3db_xy = c3db_xy_all;

if save_opt==1
    save(fullfile(save_path,[save_header,'_indiv_click_contour.mat']),'-struct','A');
end


% Plot
contour_sm_len = 5;
colorset = jet(num_freq);

fig_mf_all = figure;
set(fig_mf_all,'position',[150 150 1200 400]);
for iF=1:num_freq
    for iB=1:length(processed_files)
        click_side = raw_meas.click_side{iB};
        right_idx = find(click_side==1);
        left_idx = find(click_side==0);
        figure(fig_mf_all);
        subplot(121);
        hold on
        for iCL=1:length(left_idx)
            xy = c3db_xy_all{iB}{iF,left_idx(iCL)};
            xy_sm(:,1) = smooth(xy(:,1),contour_sm_len);
            xy_sm(:,2) = smooth(xy(:,2),contour_sm_len);
            xy_sm(isnan(xy(:,1)),:) = NaN;
            plot(xy_sm(:,1),xy_sm(:,2),'linewidth',0.5,'color',colorset(iF,:));
            clear xy_sm
        end
        figure(fig_mf_all);
        subplot(122);
        hold on
        for iCR=1:length(right_idx)
            xy = c3db_xy_all{iB}{iF,right_idx(iCR)};
            xy_sm(:,1) = smooth(xy(:,1),contour_sm_len);
            xy_sm(:,2) = smooth(xy(:,2),contour_sm_len);
            xy_sm(isnan(xy(:,1)),:) = NaN;
            plot(xy_sm(:,1),xy_sm(:,2),'linewidth',0.5,'color',colorset(iF,:));
            clear xy_sm
        end
    end
end
subplot(121)
axis equal
colormap(jet(num_freq))
colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
         'TickLabels',{num2str(freq_wanted'/1e3)},'location','southoutside');
axis([map_plot_xlim,map_plot_ylim])
grid
subplot(122)
axis equal
colormap(jet(num_freq))
colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
         'TickLabels',{num2str(freq_wanted'/1e3)},'location','southoutside');
axis([map_plot_xlim,map_plot_ylim])
grid

if save_opt==1
    saveSameSize(fig_mf_all,'file',fullfile(save_path,[save_header,'indiv_click_contour.png']),...
        'format','png','renderer','painters');
end
warning on
