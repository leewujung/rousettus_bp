% 2015 12 03  Merge data from all clicks
% 2015 12 12  Use 35 kHz rotation for other frequencies
% 2015 12 22  Adopt the merging click code for finding multi-freq center of
%             the beampatterns0

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
        [click_side,raw,rot_max,rot_elps_ctr,rot_elps_ctr_tilt,fig_elp] = shift_rotate_bp(az(ch_include_idx),el(ch_include_idx),call_dB(ch_include_idx),'eckert4',0);

        az_max_all{iB}(iC_save,:) = rot_max.az;
        el_max_all{iB}(iC_save,:) = rot_max.el;
        
        az_ectr_all{iB}(iC_save,:) = rot_elps_ctr.az;
        el_ectr_all{iB}(iC_save,:) = rot_elps_ctr.el;
        x_ectr_all{iB}(iC_save,:) = rot_elps_ctr.x;
        y_ectr_all{iB}(iC_save,:) = rot_elps_ctr.y;
        
        az_ectr_tilt_all{iB}(iC_save,:) = rot_elps_ctr_tilt.az;
        el_ectr_tilt_all{iB}(iC_save,:) = rot_elps_ctr_tilt.el;

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
            [C,~] = contour(rot_elps_ctr_tilt.xq,rot_elps_ctr_tilt.yq,raw.vq_norm,0:-3:-39,'fill','on');
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



%% Plot contours of the same frequency across all clicks together
[xqlim,yqlim] = mfwdtran(mstruct,[-90 90 0 0],[0 0 -180 180]);
xqlim = xqlim(3:4);
yqlim = yqlim(1:2);
colorset = jet(num_freq);

[freq_wanted_sort,freq_sort_idx] = sort(freq_wanted);
sm_len = 5;

fig_mf_all = figure;
set(fig_mf_all,'position',[150 150 1200 400]);
for iB=1:length(bat_proc_file_all)
    right_idx = find(click_side_single_all{iB}==1);
    left_idx = find(click_side_single_all{iB}==0);
    for iF=1:length(freq_wanted_sort)
        iF_save = freq_sort_idx(iF);
        figure(fig_mf_all);
        subplot(121);
        hold on
        for iCL=1:length(left_idx)
            xy = c3db_xy_all{iB}{iF_save,left_idx(iCL)};
            xy_sm(:,1) = smooth(xy(:,1),sm_len);
            xy_sm(:,2) = smooth(xy(:,2),sm_len);
            xy_sm(isnan(xy(:,1)),:) = NaN;
            plot(xy_sm(:,1),xy_sm(:,2),'linewidth',2,'color',colorset(iF,:));
            clear xy_sm
        end
        figure(fig_mf_all);
        subplot(122);
        hold on
        for iCR=1:length(right_idx)
            xy = c3db_xy_all{iB}{iF_save,right_idx(iCR)};
            xy_sm(:,1) = smooth(xy(:,1),sm_len);
            xy_sm(:,2) = smooth(xy(:,2),sm_len);
            xy_sm(isnan(xy(:,1)),:) = NaN;
            plot(xy_sm(:,1),xy_sm(:,2),'linewidth',2,'color',colorset(iF,:));
            clear xy_sm
        end
    end
end
subplot(121)
axis equal
colormap(jet(num_freq))
colorbar('Ticks',linspace(0,1,num_freq),'TickLabels',{num2str(freq_wanted_sort'/1e3)})
axis([xqlim,yqlim])
grid
subplot(122)
axis equal
colormap(jet(num_freq))
colorbar('Ticks',linspace(0,1,num_freq),'TickLabels',{num2str(freq_wanted_sort'/1e3)})
axis([xqlim,yqlim])
grid



%% Composite click analysis

% Merge all data
click_side_merge = reshape(cell2mat(click_side_all),1,[]);
click_num_total = size(click_side_merge)/length(az);
az_all = reshape(cell2mat(az_ectr_tilt_all),1,[]);  % rotate to ellipse center and compensate for rotation
el_all = reshape(cell2mat(el_ectr_tilt_all),1,[]);
click_max_x = cell2mat(elp_ctr_x_all');
click_max_y = cell2mat(elp_ctr_y_all');

[~,freq_sort_idx] = sort(freq_wanted);
click_max_x_sort = click_max_x(freq_sort_idx,:);
click_max_y_sort = click_max_y(freq_sort_idx,:);

% Combine all clicks and find contours (multi-freq)
threshold = 3;  % only use averaged results if number of measurement > threshold
binsize = 10;  % bin size in azimuth and elevation

az_right = az_all(click_side_merge==1);
idx_nan_right = isnan(az_right);  % index of NaN data point
click_num_right = size(az_right)/length(az);
az_right(idx_nan_right) = [];
el_right = el_all(click_side_merge==1);
el_right(idx_nan_right) = [];

az_left = az_all(click_side_merge==0);
idx_nan_left = isnan(az_left);  % index of NaN data point
click_num_left = size(az_left)/length(az);
az_left(idx_nan_left) = [];
el_left = el_all(click_side_merge==0);
el_left(idx_nan_left) = [];

call_dB_norm_merge = nan(num_freq,length(click_side_merge));
c3db_xy_all_right = cell(num_freq,1);
c3db_xy_all_left = cell(num_freq,1);
call_dB_right = zeros(num_freq,sum(~idx_nan_right));
call_dB_left = zeros(num_freq,sum(~idx_nan_left));
interp_avg_right(num_freq) = [];
interp_avg_left(num_freq) = [];
bin_avg_right(num_freq) = [];
bin_avg_left(num_freq) = [];

fig_mf = figure;
set(fig_mf,'position',[150 150 1200 400]);
for iF=1:num_freq
    % All data points for a single freq
    call_dB_norm_merge(iF,:) = reshape(cell2mat(call_dB_norm_all(:,iF)),1,[]);

    tmp = call_dB_norm_merge(iF,click_side_merge==1);
    tmp(idx_nan_right) = [];  % take out NaN data point
    call_dB_right(iF,:) = tmp;
    tmp = call_dB_norm_merge(iF,click_side_merge==0);
    tmp(idx_nan_left) = [];  % take out NaN data point
    call_dB_left(iF,:) = tmp;
    
    % Interpolation using all data points
%     [~,vq_norm_all_right(iF,:,:),azq_right,elq_right] = interp_bp(az_right/180*pi,el_right/180*pi,call_dB_right(iF,:),'rbf');
%     [~,vq_norm_all_left(iF,:,:),azq_left,elq_left] = interp_bp(az_left/180*pi,el_left/180*pi,call_dB_left(iF,:),'rbf');
    
    % Interpolation using averaged data
    [interp_avg_right(iF),bin_avg_right(iF)] = average_call(az_right,el_right,call_dB_right(iF,:),binsize,'eckert4',threshold);
    [interp_avg_left(iF),bin_avg_left(iF)] = average_call(az_left,el_left,call_dB_left(iF,:),binsize,'eckert4',threshold);
    
    % Get -3dB contour from averaged data
    [~,c_level_nan_left] = get_main_contour(interp_avg_left(iF).vq_norm_avg,interp_avg_left(iF).azq_avg,interp_avg_left(iF).elq_avg,-3);
    [~,c_level_nan_right] = get_main_contour(interp_avg_right(iF).vq_norm_avg,interp_avg_right(iF).azq_avg,interp_avg_right(iF).elq_avg,-3);
    c3db_xy_all_right{iF} = c_level_nan_left;
    c3db_xy_all_left{iF} = c_level_nan_right;
    
%     fig_mf;
%     subplot(121)  % -3dB contour for averaged left click
%     [C,~] = contour(interp_avg_left(iF).xq_avg,interp_avg_left(iF).yq_avg,interp_avg_left(iF).vq_norm_avg,0:-3:-39,'fill','on');
%     Cout = parse_contour_output(C);
%     c3db_xy_left = [];
%     for iT=1:length(Cout)  % in case contour break into pieces
%         if Cout(iT).Level == -3
%             c3db_xy_left = [c3db_xy_left; Cout(iT).X',Cout(iT).Y'];
%         end
%     end
%     title('Left click');
%     subplot(122)  % -3dB contour for averaged right click
%     [C,~] = contour(interp_avg_right(iF).xq_avg,interp_avg_right(iF).yq_avg,interp_avg_right(iF).vq_norm_avg,0:-3:-39,'fill','on');
%     Cout = parse_contour_output(C);
%     c3db_xy_right = [];
%     for iT=1:length(Cout)  % in case contour break into pieces
%         if Cout(iT).Level == -3
%             c3db_xy_right = [c3db_xy_right; Cout(iT).X',Cout(iT).Y'];
%         end
%     end
%     title('Right click');
%     suptitle(sprintf('%d kHz merged beampattern',freq_wanted(iF)/1e3));
%     save_fname = sprintf('freq_%dkHz.png',freq_wanted(iF)/1e3);
%     saveSameSize(fig_mf,'file',fullfile(save_path,save_fname),'format','png','renderer','painters');
%     pause(0.5);
    
end


%% Plot multi-freq contour
[xqlim,yqlim] = mfwdtran(mstruct,[-90 90 0 0],[0 0 -180 180]);
xqlim = xqlim(3:4);
yqlim = yqlim(1:2);
colorset = jet(num_freq);

[freq_wanted_sort,iF_seq] = sort(freq_wanted);

fig_all_freq = figure;
set(fig_all_freq,'position',[150 150 1200 400]);
for iF=iF_seq
    subplot(121)
    hold on
    plot(c3db_xy_all_left{iF}(:,1),c3db_xy_all_left{iF}(:,2),'linewidth',2,'color',colorset(iF,:));
    subplot(122)
    hold on
    plot(c3db_xy_all_right{iF}(:,1),c3db_xy_all_right{iF}(:,2),'linewidth',2,'color',colorset(iF,:));
end
subplot(121)
axis equal
colormap(jet(num_freq))
colorbar('Ticks',linspace(0,1,num_freq),'TickLabels',{num2str(freq_wanted_sort'/1e3)})
axis([xqlim,yqlim])
grid
subplot(122)
axis equal
colormap(jet(num_freq))
colorbar('Ticks',linspace(0,1,num_freq),'TickLabels',{num2str(freq_wanted_sort'/1e3)})
axis([xqlim,yqlim])
grid
saveSameSize(fig_all_freq,'file',fullfile(save_path,'all_freq_3dB.png'),'format','png','renderer','painters');


%% Scatter and interpolation for all frequencies
for iF=1:num_freq
    % Scatter plot showing all data points
    figure
    scatter(az_right,el_right,15,call_dB_right(iF,:),'fill')
    xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
    title(sprintf('All right click data points, %d kHz',freq_wanted(iF)/1e3));
    colorbar
    grid
    saveas(gcf,fullfile(save_path,sprintf('all_pts_freq%dkHz_right.png',freq_wanted(iF)/1e3)),'png');
    
    figure
    scatter(az_left,el_left,15,call_dB_left(iF,:),'fill')
    xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
    title(sprintf('All left click data points, %d kHz',freq_wanted(iF)/1e3));
    colorbar
    grid
    saveas(gcf,fullfile(save_path,sprintf('all_pts_freq%dkHz_left.png',freq_wanted(iF)/1e3)),'png');
    
    % No averaging, right
    figure
    axesm eckert4
    framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
    gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
    axis off
    contourfm(elq_right/pi*180,azq_right/pi*180,squeeze(vq_norm_all_right(iF,:,:)),0:-3:-30);
    caxis([-30 0])
    colormap(parula(10))
    colorbar
    title(sprintf('Merged right clicks, %d kHz, no averaging',freq_wanted(iF)/1e3));
    saveas(gcf,fullfile(save_path,sprintf('merged_no_avg_freq%dkHz_right.png',freq_wanted(iF)/1e3)),'png');
        
    % No averaging, left
    figure
    axesm eckert4
    framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
    gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
    axis off
    tmp = squeeze(vq_norm_all_left(iF,:,:));
%     contourfm(elq_left/pi*180,azq_left/pi*180,squeeze(vq_norm_all_left(iF,:,:)),0:-3:-30);
    contourfm(elq_left/pi*180,azq_left/pi*180,tmp,0:-3:-30);
    caxis([-30 0])
    colormap(parula(10))
    colorbar
    title(sprintf('Merged left clicks, %d kHz, no averaging',freq_wanted(iF)/1e3));
    saveas(gcf,fullfile(save_path,sprintf('merged_no_avg_freq%dkHz_left.png',freq_wanted(iF)/1e3)),'png');
    
    % Averaging, right
    figure
    axesm eckert4
    framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
    gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
    axis off
    contourfm(interp_avg_right(iF).elq_avg,interp_avg_right(iF).azq_avg,interp_avg_right(iF).vq_norm_avg,0:-3:-30);
    caxis([-30 0])
    colormap(parula(10))
    colorbar
    title(sprintf('Averaged right clicks, %d kHz, th=%d, bin=%ddeg',freq_wanted(iF)/1e3,threshold,binsize));
    saveas(gcf,fullfile(save_path,sprintf('merged_avg_freq%dkHz_right.png',freq_wanted(iF)/1e3)),'png');
    
    % Averaging, left
    figure
    axesm eckert4
    framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
    gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
    axis off
    contourfm(interp_avg_left(iF).elq_avg,interp_avg_left(iF).azq_avg,interp_avg_left(iF).vq_norm_avg,0:-3:-30);
    caxis([-30 0])
    colormap(parula(10))
    colorbar
    title(sprintf('Averaged left clicks, %d kHz, th=%d, bin=%ddeg',freq_wanted(iF)/1e3,threshold,binsize));
    saveas(gcf,fullfile(save_path,sprintf('merged_avg_freq%dkHz_left.png',freq_wanted(iF)/1e3)),'png');
end

