% 2015 12 03  Merge data from all clicks

clear
warning off

usrn = getenv('username');
base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
bat_proc_path = './proc_output';
bat_proc_file_all = dir(fullfile(base_path,bat_proc_path,'rousettus_20150825_36134_*'));
% bat_proc_file = 'rousettus_20150825_36134_05_mic_data_bp_proc';
plot_opt = 0;
freq_wanted = 35e3;

addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\EllipseDirectFit']);
% addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\beampattern_processing']);

az_rot_ecen_all = cell(length(bat_proc_file_all),1);
el_rot_ecen_all = az_rot_ecen_all;
az_rot_max_all = az_rot_ecen_all;
el_rot_max_all = el_rot_ecen_all;
x_rot_ecen_all = az_rot_ecen_all;
y_rot_ecen_all = az_rot_ecen_all;
call_dB_norm_all = az_rot_ecen_all;
click_side_all = az_rot_ecen_all;

for iF = 1:length(bat_proc_file_all)
%     iF = 2;
    bat_proc_file = bat_proc_file_all(iF).name;
    fprintf('File: %s\n',bat_proc_file);
    
    data = load(fullfile(base_path,bat_proc_path,bat_proc_file));
    mic_num = 1:data.mic_data.num_ch_in_file;
    good_call_idx = find(data.proc.chk_good_call);
    
    if plot_opt
        fig_all = figure;
    end
    numrow = ceil(length(good_call_idx)/2);
    
    for iC = good_call_idx'
%         iC = 14;
        iC_save = find(iC==good_call_idx);
        
        % Get call info
        [call_dB,az,el,ch_include_idx] = get_call_azel_dB_data(data,freq_wanted,iC);
        [vq,vq_norm,azq,elq] = interp_bp(az(ch_include_idx),el(ch_include_idx),call_dB(ch_include_idx),'rbf');
        call_dB_norm = call_dB-max(call_dB);
        az = az/pi*180;  % convert to [deg]
        el = el/pi*180;
        azq = azq/pi*180;
        elq = elq/pi*180;
        
        % Set map projection
        mstruct = defaultm('ortho');
        mstruct = defaultm(mstruct);

        % Rotate measurements to use max position as origin
        [mm,mmidx] = max(vq_norm(:));
        origin = [elq(mmidx),azq(mmidx)];  % [Lat Lon]
        [el_rot_max,az_rot_max] = rotatem(el,az,origin,'forward','degrees');
        [elq_rot_max,azq_rot_max] = rotatem(elq,azq,origin,'forward','degrees');
        
        % Determine right/left click
        if azq(mmidx)<0
            click_side = zeros(length(az),1);  % left: click_side = 0
            click_side_t = 'left';
        else
            click_side = ones(length(az),1);  % right: click_side = 1
            click_side_t = 'right';
        end

        % Project lat-lon to map projection distance
        [xq,yq] = mfwdtran(mstruct,elq,azq);  % rotation
        [xq_rot_max,yq_rot_max] = mfwdtran(mstruct,elq_rot_max,azq_rot_max);  % map projection

        % Fit ellipse to projected beampattern
        fig_elp = figure;
        subplot(121)
        [C,~] = contour(xq_rot_max,yq_rot_max,vq_norm,0:-3:-39,'fill','on');
        Cout = parse_contour_output(C);
        c3db_xy = [];
        for iT=1:length(Cout)  % in case contour break into pieces
            if Cout(iT).Level == -3
                c3db_xy = [c3db_xy; Cout(iT).X',Cout(iT).Y'];
            end
        end
        axis([-pi pi -pi/2 pi/2]);
        axis equal
        
        subplot(122)
        contour(xq,yq,vq_norm,0:-3:-39,'fill','on');
        axis([-pi pi -pi/2 pi/2]);
        axis equal
        title('Not rotated');
        
        A = EllipseDirectFit(c3db_xy);  % fit ellipse (direct fit)
        E = get_ellipse_param(A);       % get ellipse parameters
        
        % Plot ellipse
        xmin = min(c3db_xy(:,1))-range(c3db_xy(:,1)*0.5);
        xmax = max(c3db_xy(:,1))+range(c3db_xy(:,1)*0.5);
        ymin = min(c3db_xy(:,2))-range(c3db_xy(:,1)*0.5);
        ymax = max(c3db_xy(:,2))+range(c3db_xy(:,1)*0.5);
        figure(fig_elp)
        subplot(121)
        hold on
        fit_df = ezplot(E.eqt,[xmin,xmax,ymin,ymax]);
        set(fit_df,'linecolor','b','linewidth',2);
        hold off
        title(sprintf('Call #%02d, ar=%2.2f, theta=%2.2fdeg, side=%s',iC,E.ar,E.theta/pi*180,click_side_t));
        close
        
        % Rotate measurements to use ellipse center as origin
        [lat_ecen,lon_ecen] = minvtran(mstruct,E.x0,E.y0);
        [el_rot_ecen,az_rot_ecen] = rotatem(el,az,[lat_ecen,lon_ecen],'forward','degrees');  % rotation
        [elq_rot_ecen,azq_rot_ecen] = rotatem(elq,azq,[lat_ecen,lon_ecen],'forward','degrees');
        [x_rot_ecen,y_rot_ecen] = mfwdtran(mstruct,el_rot_ecen,az_rot_ecen);  % map projection
        [xq_rot_ecen,yq_rot_ecen] = mfwdtran(mstruct,elq_rot_ecen,azq_rot_ecen);
        
        az_rot_max_all{iF}(iC_save,:) = az_rot_max;
        el_rot_max_all{iF}(iC_save,:) = el_rot_max;
        az_rot_ecen_all{iF}(iC_save,:) = az_rot_ecen;
        el_rot_ecen_all{iF}(iC_save,:) = el_rot_ecen;
        x_rot_ecen_all{iF}(iC_save,:) = x_rot_ecen;
        y_rot_ecen_all{iF}(iC_save,:) = y_rot_ecen;
        call_dB_norm_all{iF}(iC_save,:) = call_dB_norm;
        click_side_all{iF}(iC_save,:) = click_side;
        
        
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


% Merge all data
click_side_merge = reshape(cell2mat(click_side_all),1,[]);
% az_all = reshape(cell2mat(az_rot_ecen_all),1,[]);
% el_all = reshape(cell2mat(el_rot_ecen_all),1,[]);
az_all = reshape(cell2mat(az_rot_max_all),1,[]);
el_all = reshape(cell2mat(el_rot_max_all),1,[]);
x_all = reshape(cell2mat(x_rot_ecen_all),1,[]);
y_all = reshape(cell2mat(y_rot_ecen_all),1,[]);
call_dB_norm_merge = reshape(cell2mat(call_dB_norm_all),1,[]);


% All right clicks
az_right = az_all(click_side_merge==1);
el_right = el_all(click_side_merge==1);
call_dB_right = call_dB_norm_merge(click_side_merge==1);
az_left = az_all(click_side_merge==0);
el_left = el_all(click_side_merge==0);
call_dB_left = call_dB_norm_merge(click_side_merge==0);
threshold = 3;

[~,vq_norm_all_right,azq_right,elq_right] = interp_bp(az_right/180*pi,el_right/180*pi,call_dB_right,'rbf');
[~,vq_norm_all_left,azq_left,elq_left] = interp_bp(az_left/180*pi,el_left/180*pi,call_dB_left,'rbf');

binsize = 10;
[interp_avg_right,bin_avg_right] = average_call(az_right,el_right,call_dB_right,binsize,'eckert4',threshold);
[interp_avg_left,bin_avg_left] = average_call(az_left,el_left,call_dB_left,binsize,'eckert4',threshold);


%% Plot
% Scatter plot showing all data points
figure
scatter(az_right,el_right,15,call_dB_right,'fill')
xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
title('All right click data points');
colorbar
grid

figure
scatter(az_left,el_left,15,call_dB_left,'fill')
xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
title('All left click data points');
colorbar
grid


% No averaging, right
figure
axesm eckert4
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
axis off
contourfm(elq_right/pi*180,azq_right/pi*180,vq_norm_all_right,0:-3:-30);
caxis([-30 0])
colormap(parula(10))
colorbar
title(sprintf('Averaged right clicks, th=%d, bin=%ddeg',threshold,binsize));

% No averaging, left
figure
axesm eckert4
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
axis off
contourfm(elq_left/pi*180,azq_left/pi*180,vq_norm_all_left,0:-3:-30);
caxis([-30 0])
colormap(parula(10))
colorbar
title(sprintf('Averaged left clicks, th=%d, bin=%ddeg',threshold,binsize));

% Averaging, right
figure
axesm eckert4
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
axis off
contourfm(interp_avg_right.elq_avg,interp_avg_right.azq_avg,interp_avg_right.vq_norm_avg,0:-3:-30);
caxis([-30 0])
colormap(parula(10))
colorbar
title(sprintf('Averaged right clicks, th=%d, bin=%ddeg',threshold,binsize));

% Averaging, left
figure
axesm eckert4
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
axis off
contourfm(interp_avg_left.elq_avg,interp_avg_left.azq_avg,interp_avg_left.vq_norm_avg,0:-3:-30);
caxis([-30 0])
colormap(parula(10))
colorbar
title(sprintf('Averaged left clicks, th=%d, bin=%ddeg',threshold,binsize));




