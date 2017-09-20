% 2015 12 23  Test contour plot

clear
load('C:\Users\Wu-Jung Lee\Dropbox\0_ANALYSIS\bp_processing\20151216_merge_clicks_multifreq\all_click_data.mat');

% Merge all data
click_side_merge = reshape(cell2mat(click_side_all),1,[]);
click_num_total = size(click_side_merge)/length(az);
az_all = reshape(cell2mat(az_ecen_tilt_all),1,[]);  % rotate to ellipse center and compensate for rotation
el_all = reshape(cell2mat(el_ecen_tilt_all),1,[]);

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

freq_wanted_sort = sort(freq_wanted);
[xqlim,yqlim] = mfwdtran(mstruct,[-90 90 0 0],[0 0 -180 180]);
xqlim = xqlim(3:4);
yqlim = yqlim(1:2);
colorset = jet(num_freq);

fig_mf = figure;
set(fig_mf,'position',[150 150 1200 400]);

for iF=1:length(freq_wanted)
    iF_color = find(freq_wanted(iF)==freq_wanted_sort);
    
    % All data points for a single freq
    call_dB_norm_merge(iF,:) = reshape(cell2mat(call_dB_norm_all(:,iF)),1,[]);
    
    tmp = call_dB_norm_merge(iF,click_side_merge==1);
    tmp(idx_right) = [];
    call_dB_right(iF,:) = tmp;
    tmp = call_dB_norm_merge(iF,click_side_merge==0);
    tmp(idx_left) = [];
    call_dB_left(iF,:) = tmp;
    
    % Interpolation using averaged data
    [interp_avg_right(iF),bin_avg_right(iF)] = average_call(az_right,el_right,call_dB_right(iF,:),binsize,'eckert4',threshold);
    [interp_avg_left(iF),bin_avg_left(iF)] = average_call(az_left,el_left,call_dB_left(iF,:),binsize,'eckert4',threshold);
    
    % Get -3dB contour from averaged data
    C_main_left = get_main_contour(interp_avg_left(iF).vq_norm_avg,interp_avg_left(iF).azq_avg(1,:),interp_avg_left(iF).elq_avg(:,1),-3);
    Cout_left = parse_contour_output(C_main_left);
    [c3db_left_xy(:,1),c3db_left_xy(:,2)] = mfwdtran(mstruct,Cout_left.Y,Cout_left.X);
    C_main_right = get_main_contour(interp_avg_right(iF).vq_norm_avg,interp_avg_right(iF).azq_avg(1,:),interp_avg_right(iF).elq_avg(:,1),-3);
    Cout_right = parse_contour_output(C_main_right);
    [c3db_right_xy(:,1),c3db_right_xy(:,2)] = mfwdtran(mstruct,Cout_right.Y,Cout_right.X);
    
    sm_len = 5;
    
    figure(fig_mf);
    subplot(121)
    hold on
    plot(smooth(c3db_left_xy(:,1),sm_len),smooth(c3db_left_xy(:,2),sm_len),'linewidth',2,'color',colorset(iF_color,:));
    subplot(122)
    hold on
    plot(smooth(c3db_right_xy(:,1),sm_len),smooth(c3db_right_xy(:,2),sm_len),'linewidth',2,'color',colorset(iF_color,:));
    
    clear c3db_*_xy
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
