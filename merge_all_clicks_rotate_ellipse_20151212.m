% 2015 12 03  Merge data from all clicks
% 2015 12 12  Use 35 kHz rotation for other frequencies

clear
warning off

usrn = getenv('username');
base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
save_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20151216_merge_clicks_multifreq'];
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
az_ecen_all = cell(length(bat_proc_file_all),num_freq);
el_ecen_all = az_ecen_all;
az_max_all = az_ecen_all;
el_max_all = az_ecen_all;
az_ecen_tilt_all = az_ecen_all;
el_ecen_tilt_all = az_ecen_all;
x_ecen_all = az_ecen_all;
y_ecen_all = az_ecen_all;
call_dB_norm_all = az_ecen_all;
click_side_all = az_ecen_all;

for iB = 1:length(bat_proc_file_all)
%     iF = 2;
    bat_proc_file = bat_proc_file_all(iB).name;
    fprintf('File: %s\n',bat_proc_file);
    
    data = load(fullfile(base_path,bat_proc_path,bat_proc_file));
    mic_num = 1:data.mic_data.num_ch_in_file;
    good_call_idx = find(data.proc.chk_good_call);
    
    if plot_opt
        fig_all = figure;
    end
    numrow = ceil(length(good_call_idx)/2);
    
    for iC = good_call_idx'
%         iC = 11;
        iC_save = find(iC==good_call_idx);
        
        % Get call info
        for iF=1:num_freq
            [call_dB(iF,:),az,el,ch_include_idx] = get_call_azel_dB_data(data,freq_wanted(iF),iC);
        end
        [~,vq_norm,azq,elq] = interp_bp(az(ch_include_idx),el(ch_include_idx),call_dB(1,ch_include_idx),'rbf');  % use the first frequency data for finding ellipse center
        call_dB_norm = bsxfun(@minus,call_dB,max(call_dB,[],2));
        az = az/pi*180;  % convert to [deg]
        el = el/pi*180;
        azq = azq/pi*180;
        elq = elq/pi*180;
        
        % Rotate measurements to use max position as origin
        [mm,mmidx] = max(vq_norm(:));
        origin_max = [elq(mmidx),azq(mmidx)];  % [Lat Lon]
        [el_max,az_max] = rotatem(el,az,origin_max,'forward','degrees');
        [elq_max,azq_max] = rotatem(elq,azq,origin_max,'forward','degrees');
        
        % Determine right/left click
        if azq(mmidx)<0
            click_side = zeros(length(az),1);  % left: click_side = 0
            click_side_t = 'left';
        else
            click_side = ones(length(az),1);  % right: click_side = 1
            click_side_t = 'right';
        end

        % Project lat-lon to map projection distance
        [xq,yq] = mfwdtran(mstruct,elq,azq);
        [x_max,y_max] = mfwdtran(mstruct,el_max,az_max);  % map projection
        [xq_max,yq_max] = mfwdtran(mstruct,elq_max,azq_max);  % map projection

        % Plot raw data
        fig_elp = figure;
        set(fig_elp,'position',[150 150 1200 400]);
        subplot(141)
        contour(xq,yq,vq_norm,0:-3:-39,'fill','on');
        axis equal; axis([-1.1 1.1 -1.1 1.1]); grid on
        title('Raw data');
        
        % Plot measurements rotated to max point and fit ellipse
        E_max = plot_bp_fit_ellipse(subplot(142),xq_max,yq_max,vq_norm);
        title('Best-fitting ellipse');

        % Iteratively rotate/shift until the center of best-fitting ellipse
        % is at [0,0] in X-Y plane
        fig_it = figure;
        aax = gca;
        n = 1;
        E1 = E_max;
        x1 = x_max;       y1 = y_max;
        el1 = el_max;     az1 = az_max;
        xq1 = xq_max;     yq1 = yq_max;
        elq1 = elq_max;   azq1 = azq_max;
        while sqrt(E1.x0^2+E1.y0^2)>0.01
            E1 = plot_bp_fit_ellipse(aax,xq1,yq1,vq_norm);
            [origin1(1),origin1(2)] = minvtran(mstruct,E1.x0,E1.y0);  % transform ellipse center x-y pos to lat-lon
            [elq1,azq1] = rotatem(elq1,azq1,origin1,'forward','degrees');  % rotation
            [el1,az1] = rotatem(el1,az1,origin1,'forward','degrees');  % rotation
            [xq1,yq1] = mfwdtran(mstruct,elq1,azq1);  % map projection
            [x1,y1] = mfwdtran(mstruct,el1,az1);  % map projection
            n = n+1;
%             pause(0.5)
        end
        close(fig_it)
        fprintf('Rotate %d times\n',n);
        x_ecen = x1;      y_ecen = y1;
        el_ecen = el1;    az_ecen = az1;
        xq_ecen = xq1;    yq_ecen = yq1;
        elq_ecen = elq1;  azq_ecen = azq1;
        E_ecen = E1;
        figure(fig_elp)
        E_ecen2 = plot_bp_fit_ellipse(subplot(143),xq_ecen,yq_ecen,vq_norm);
        title(sprintf('Tilt %2.2fdeg',E_ecen.theta/pi*180));
        text(-1,-1,sprintf('a0=%2.2f, b0=%2.2f',E_ecen.a0,E_ecen.b0));
        
        % Compensate for ellipse tilt
        if E_ecen.theta<0 && E_ecen.coef.a<E_ecen.coef.c && E_ecen.a0>E_ecen.b0  % horizontal ellipse
            elp_theta = -(pi/2+E_ecen.theta);
        elseif E_ecen.theta>0 && E_ecen.coef.a<E_ecen.coef.c && E_ecen.a0>E_ecen.b0  % horizontal ellipse
            elp_theta = pi/2-E_ecen.theta;
        else
            elp_theta = -E_ecen.theta;
        end
        
        xq_ecen_tilt = xq_ecen*cos(elp_theta) - yq_ecen*sin(elp_theta);
        yq_ecen_tilt = xq_ecen*sin(elp_theta) + yq_ecen*cos(elp_theta);
        x_ecen_tilt = x_ecen*cos(elp_theta) - y_ecen*sin(elp_theta);
        y_ecen_tilt = x_ecen*sin(elp_theta) + y_ecen*cos(elp_theta);
        
        [xxx,yyy] = mfwdtran(mstruct,[0,0,-90,90],[-180,180,0,0]);
        xxx = xxx(1:2);
        yyy = yyy(3:4);
        idx_right = x_ecen_tilt>xxx(1) & x_ecen_tilt<xxx(2) & y_ecen_tilt>yyy(1) & y_ecen_tilt<yyy(2);
        el_ecen_tilt = nan(length(x_ecen_tilt),1);
        az_ecen_tilt = nan(length(y_ecen_tilt),1);
        [el_ecen_tilt(idx_right),az_ecen_tilt(idx_right)] = minvtran(mstruct,x_ecen_tilt(idx_right),y_ecen_tilt(idx_right));

        figure(fig_elp)
        E_ecen_tilt = plot_bp_fit_ellipse(subplot(144),xq_ecen_tilt,yq_ecen_tilt,vq_norm);  % check tilt compensation results
        title('Tilt compensated')

        suptitle(sprintf('Call #%02d, ar=%2.2f, side=%s',iC,E_ecen.ar,click_side_t));
%         save_fname = sprintf('file%d_call%d.png',iF,iC);
%         saveSameSize(fig_elp,'file',fullfile(save_path,save_fname),'format','png','renderer','painters');
%         pause(0.5)
        close(fig_elp)
        
        az_max_all{iB}(iC_save,:) = az_max;
        el_max_all{iB}(iC_save,:) = el_max;
        x_max_all{iB}(iC_save,:) = x_max;
        y_max_all{iB}(iC_save,:) = y_max;
        
        az_ecen_all{iB}(iC_save,:) = az_ecen;
        el_ecen_all{iB}(iC_save,:) = el_ecen;
        x_ecen_all{iB}(iC_save,:) = x_ecen;
        y_ecen_all{iB}(iC_save,:) = y_ecen;
        
        az_ecen_tilt_all{iB}(iC_save,:) = az_ecen_tilt;
        el_ecen_tilt_all{iB}(iC_save,:) = el_ecen_tilt;
        x_ecen_tilt_all{iB}(iC_save,:) = x_ecen_tilt;
        y_ecen_tilt_all{iB}(iC_save,:) = y_ecen_tilt;

        for iF=1:num_freq
            call_dB_norm_all{iB,iF}(iC_save,:) = call_dB_norm(iF,:);
        end
        click_side_all{iB}(iC_save,:) = click_side;
        
        
        % Plot rotated bp using ellipse center
        if plot_opt
            figure(fig_all)
            subplot(2,numrow,find(iC==good_call_idx));
            [C,~] = contour(xq_rot_ecen,yq_rot_ecen,vq_norm,0:-3:-39,'fill','on');
            title(sprintf('Call #%02d',iC));
        end
    end
    
    if plot_opt
        figure(fig_all)
        suptitle(bat_proc_file);
    end
    
end
warning on

save(fullfile(save_path,'all_click_data.mat'));


% Merge all data
click_side_merge = reshape(cell2mat(click_side_all),1,[]);
click_num_total = size(click_side_merge)/length(az);
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
    [C,h] = contour(interp_avg_left(iF).xq_avg,interp_avg_left(iF).yq_avg,interp_avg_left(iF).vq_norm_avg,0:-3:-39,'fill','on');
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
%     save_fname = sprintf('freq_%dkHz.png',freq_wanted(iF)/1e3);
%     saveSameSize(fig_mf,'file',fullfile(save_path,save_fname),'format','png','renderer','painters');
%     pause(0.5);
    
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
