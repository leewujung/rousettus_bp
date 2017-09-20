% 2015 12 16  Center of 3dB beampattern across frequency
% 2015 12 22  Center of shifted beampattern

clear
warning off

usrn = getenv('username');
base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
save_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20151216_multifreq_trace'];
if ~exist(save_path,'dir')
    mkdir(save_path);
end
bat_proc_path = './proc_output_rousettus_checked';
bat_proc_file_all = dir(fullfile(base_path,bat_proc_path,'rousettus_20150825_36134_*'));
plot_opt = 0;
freq_wanted = [20:5:50]*1e3;
num_freq = length(freq_wanted);

addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\EllipseDirectFit']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);


% Set map projection
mstruct = defaultm('eckert4');
mstruct = defaultm(mstruct);
        
% Process all files
click_side_all = cell(length(bat_proc_file_all),1);
call_max_loc_az_all = cell(length(bat_proc_file_all),1);
call_max_loc_el_all = cell(length(bat_proc_file_all),1);
call_max_loc_x_all = cell(length(bat_proc_file_all),1);
call_max_loc_y_all = cell(length(bat_proc_file_all),1);

for iB = 1:5%length(bat_proc_file_all)
%     iF = 2;
    bat_proc_file = bat_proc_file_all(iB).name;
    fprintf('File: %s\n',bat_proc_file);
    
    data = load(fullfile(base_path,bat_proc_path,bat_proc_file));
    mic_num = 1:data.mic_data.num_ch_in_file;
    good_call_idx = find(data.proc.chk_good_call);
    
    call_dB = nan(num_freq,data.mic_data.num_ch_in_file);
    for iC = good_call_idx'
%         iC = 11;
        iC_save = find(iC==good_call_idx);
        
        % Get call info
        call_max_loc_x = nan(num_freq,1);
        call_max_loc_y = nan(num_freq,1);
        fig_elp = figure;
        for iF=1:num_freq
            [call_dB(iF,:),az,el,ch_include_idx] = get_call_azel_dB_data(data,freq_wanted(iF),iC);
            [~,vq_norm,azq,elq] = interp_bp(az(ch_include_idx),el(ch_include_idx),call_dB(iF,ch_include_idx),'rbf');  % use the first frequency data for finding ellipse center
            [xq,yq] = mfwdtran(mstruct,elq/pi*180,azq/pi*180);
            E_max = plot_bp_fit_ellipse(gca,xq,yq,vq_norm);
%             pause(0.5)
            call_max_loc_x(iF,:) = E_max.x0;
            call_max_loc_y(iF,:) = E_max.y0;
        end
        close(fig_elp);
        call_max_loc_x_all{iB}(:,iC_save) = call_max_loc_x;
        call_max_loc_y_all{iB}(:,iC_save) = call_max_loc_y;

        % Determine right/left click
        [~,vq_norm,azq,elq] = interp_bp(az(ch_include_idx),el(ch_include_idx),call_dB(4,ch_include_idx),'rbf');  % use the first frequency data for finding ellipse center
        call_dB_norm = bsxfun(@minus,call_dB,max(call_dB,[],2));
        az = az/pi*180;  % convert to [deg]
        el = el/pi*180;
        [mm,mmidx] = max(vq_norm(:));
        if azq(mmidx)<0
            click_side = 0;  % left: click_side = 0
            click_side_t = 'left';
        else
            click_side = 1;  % right: click_side = 1
            click_side_t = 'right';
        end

        click_side_all{iB}(iC_save,:) = click_side;
    end
        
end
warning on

save(fullfile(save_path,'test_data.mat'));

% Merge all data
click_side_merge = reshape(cell2mat(click_side_all),1,[]);
click_num_total = size(click_side_merge)/length(az);
click_max_x = cell2mat(call_max_loc_x_all');
click_max_y = cell2mat(call_max_loc_y_all');

% az_all = reshape(cell2mat(az_max_all),1,[]);  % rotate to max mic
% el_all = reshape(cell2mat(el_max_all),1,[]);
% az_all = reshape(cell2mat(az_ecen_all),1,[]);  % rotate to ellipse center
% el_all = reshape(cell2mat(el_ecen_all),1,[]);
az_all = reshape(cell2mat(az_ecen_tilt_all),1,[]);  % rotate to ellipse center and compensate for rotation
el_all = reshape(cell2mat(el_ecen_tilt_all),1,[]);
% for iF=1:num_freq
%     call_dB_norm_merge(iF,:) = reshape(cell2mat(call_dB_norm_all),1,[]);
% end

% All right clicks
threshold = 3;
binsize = 10;

az_right = az_all(click_side_merge==1);
idx_right = isnan(az_right);
click_num_right = size(az_right)/length(az);
az_right(idx_right) = [];
el_right = el_all(click_side_merge==1);
el_right(idx_right) = [];

az_left = az_all(click_side_merge==0);
idx_left = isnan(az_left);
click_num_left = size(az_left)/length(az);
az_left(idx_left) = [];
el_left = el_all(click_side_merge==0);
el_left(idx_left) = [];

call_dB_norm_merge = nan(num_freq,length(click_side_merge));
c3db_xy_all_right = cell(num_freq,1);
c3db_xy_all_left = cell(num_freq,1);
fig_mf = figure;
set(fig_mf,'position',[150 150 1200 400]);
for iF=1:num_freq
    % All data points for a single freq
    call_dB_norm_merge(iF,:) = reshape(cell2mat(call_dB_norm_all(:,iF)),1,[]);

    tmp = call_dB_norm_merge(iF,click_side_merge==1);
    tmp(idx_right) = [];
    call_dB_right(iF,:) = tmp;
    tmp = call_dB_norm_merge(iF,click_side_merge==0);
    tmp(idx_left) = [];
    call_dB_left(iF,:) = tmp;
    
    % Interpolation using all data points
%     [~,vq_norm_all_right(iF,:,:),azq_right,elq_right] = interp_bp(az_right/180*pi,el_right/180*pi,call_dB_right(iF,:),'rbf');
%     [~,vq_norm_all_left(iF,:,:),azq_left,elq_left] = interp_bp(az_left/180*pi,el_left/180*pi,call_dB_left(iF,:),'rbf');
    
    % Interpolation using averaged data
    [interp_avg_right(iF),bin_avg_right(iF)] = average_call(az_right,el_right,call_dB_right(iF,:),binsize,'eckert4',threshold);
    [interp_avg_left(iF),bin_avg_left(iF)] = average_call(az_left,el_left,call_dB_left(iF,:),binsize,'eckert4',threshold);
    
    % Get -3dB contour from averaged data
    fig_mf;
    subplot(121)  % -3dB contour for averaged left click
    [C,~] = contour(interp_avg_left(iF).xq_avg,interp_avg_left(iF).yq_avg,interp_avg_left(iF).vq_norm_avg,0:-3:-39,'fill','on');
    Cout = parse_contour_output(C);
    c3db_xy_left = [];
    for iT=1:length(Cout)  % in case contour break into pieces
        if Cout(iT).Level == -3
            c3db_xy_left = [c3db_xy_left; Cout(iT).X',Cout(iT).Y'];
        end
    end
    title('Left click');
    subplot(122)  % -3dB contour for averaged right click
    [C,~] = contour(interp_avg_right(iF).xq_avg,interp_avg_right(iF).yq_avg,interp_avg_right(iF).vq_norm_avg,0:-3:-39,'fill','on');
    Cout = parse_contour_output(C);
    c3db_xy_right = [];
    for iT=1:length(Cout)  % in case contour break into pieces
        if Cout(iT).Level == -3
            c3db_xy_right = [c3db_xy_right; Cout(iT).X',Cout(iT).Y'];
        end
    end
    title('Right click');
    suptitle(sprintf('%d kHz merged beampattern',freq_wanted(iF)/1e3));
    save_fname = sprintf('freq_%dkHz.png',freq_wanted(iF)/1e3);
    saveSameSize(fig_mf,'file',fullfile(save_path,save_fname),'format','png','renderer','painters');
    pause(0.5);
    
    c3db_xy_all_right{iF} = c3db_xy_right;
    c3db_xy_all_left{iF} = c3db_xy_left;
    
end


% Plot
% Combined contour for all frequencies
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



% Contour plot for each frequency
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



% %% Plot
% % Scatter plot showing all data points
% figure
% scatter(az_right,el_right,15,call_dB_right,'fill')
% xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
% title('All right click data points');
% colorbar
% grid
% 
% figure
% scatter(az_left,el_left,15,call_dB_left,'fill')
% xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
% title('All left click data points');
% colorbar
% grid
% 
% 
% % No averaging, right
% figure
% axesm eckert4
% framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
% gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
% axis off
% contourfm(elq_right/pi*180,azq_right/pi*180,vq_norm_all_right,0:-3:-30);
% caxis([-30 0])
% colormap(parula(10))
% colorbar
% title(sprintf('Averaged right clicks, th=%d, bin=%ddeg',threshold,binsize));
% 
% % No averaging, left
% figure
% axesm eckert4
% framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
% gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
% axis off
% contourfm(elq_left/pi*180,azq_left/pi*180,vq_norm_all_left,0:-3:-30);
% caxis([-30 0])
% colormap(parula(10))
% colorbar
% title(sprintf('Averaged left clicks, th=%d, bin=%ddeg',threshold,binsize));
% 
% % Averaging, right
% figure
% axesm eckert4
% framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
% gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
% axis off
% contourfm(interp_avg_right.elq_avg,interp_avg_right.azq_avg,interp_avg_right.vq_norm_avg,0:-3:-30);
% caxis([-30 0])
% colormap(parula(10))
% colorbar
% title(sprintf('Averaged right clicks, th=%d, bin=%ddeg',threshold,binsize));
% 
% % Averaging, left
% figure
% axesm eckert4
% framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
% gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
% axis off
% contourfm(interp_avg_left.elq_avg,interp_avg_left.azq_avg,interp_avg_left.vq_norm_avg,0:-3:-30);
% caxis([-30 0])
% colormap(parula(10))
% colorbar
% title(sprintf('Averaged left clicks, th=%d, bin=%ddeg',threshold,binsize));
% 
% % Averaging, right
% figure
% axesm eckert4
% framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
% gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
% axis off
% contourfm(interp_avg_50k_right.elq_avg,interp_avg_50k_right.azq_avg,interp_avg_50k_right.vq_norm_avg,0:-3:-30);
% caxis([-30 0])
% colormap(parula(10))
% colorbar
% title(sprintf('Averaged right clicks, th=%d, bin=%ddeg',threshold,binsize));
% 
% % Averaging, left
% figure
% axesm eckert4
% framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
% gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
% axis off
% contourfm(interp_avg_50k_left.elq_avg,interp_avg_50k_left.azq_avg,interp_avg_50k_left.vq_norm_avg,0:-3:-30);
% caxis([-30 0])
% colormap(parula(10))
% colorbar
% title(sprintf('Averaged left clicks, th=%d, bin=%ddeg',threshold,binsize));
% 
% % Averaging, right
% figure
% axesm eckert4
% framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
% gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
% axis off
% contourfm(interp_avg_20k_right.elq_avg,interp_avg_50k_right.azq_avg,interp_avg_20k_right.vq_norm_avg,0:-3:-30);
% caxis([-30 0])
% colormap(parula(10))
% colorbar
% title(sprintf('Averaged right clicks, th=%d, bin=%ddeg',threshold,binsize));
% 
% % Averaging, left
% figure
% axesm eckert4
% framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
% gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
% axis off
% contourfm(interp_avg_20k_left.elq_avg,interp_avg_50k_left.azq_avg,interp_avg_20k_left.vq_norm_avg,0:-3:-30);
% caxis([-30 0])
% colormap(parula(10))
% colorbar
% title(sprintf('Averaged left clicks, th=%d, bin=%ddeg',threshold,binsize));
% 
% 
% 
