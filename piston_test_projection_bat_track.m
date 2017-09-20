% 2015 11 02  Use piston model to test projection onto az-el plane
% 2015 11 04  Make beam axis independent of bat head aim
% 2015 11 05/11  Use real mic location and bat tracks

clear
usrn = getenv('username');
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\beampattern_processing']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\fitellipse']);

% Bat data path
save_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_model'];
base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
bat_proc_dir = './proc_output';
bat_proc_file = dir(fullfile(base_path,bat_proc_dir,'rousettus_20150825_*.mat'));
freq = 35e3;
goodcall_angle_range = [-60 60];
mic_area_frac = 0.5;

% Load bat track and call location
A = load(fullfile(base_path,bat_proc_dir,bat_proc_file(2).name));
track = nan(length(bat_proc_file),size(A.track.track_interp,1),size(A.track.track_interp,2));
call_loc = cell(length(bat_proc_file),1);
call_loc_idx_on_track = cell(length(bat_proc_file),1);
call_max_ch_idx = cell(length(bat_proc_file),1);
call_max_azel = cell(length(bat_proc_file),1);
call_all_azel = cell(length(bat_proc_file),1);
for iB = 1:length(bat_proc_file)
    A = load(fullfile(base_path,bat_proc_dir,bat_proc_file(iB).name));
    track(iB,:,:) = A.track.track_interp;
%     good_call_idx = find(A.proc.chk_bad_call==0);
    good_call_idx = 1:length(A.proc.chk_bad_call);
    call_loc{iB} = A.proc.bat_loc_at_call(good_call_idx,:);
    call_loc_idx_on_track{iB} = A.proc.call_loc_on_track_interp(good_call_idx);
    for iC = 1:length(good_call_idx) % go through all calls to find beam aim
        [~,fidx] = min(abs(freq-A.proc.call_freq_vec{iC}));
        call_dB = squeeze(A.proc.call_psd_dB_comp_re20uPa_withbp{iC}(fidx,:));
        [~,azel_idx] = max(call_dB);
        call_max_ch_idx{iB}(iC) = azel_idx;  % maximum channel
        call_max_azel{iB}(iC,:) = A.proc.mic_to_bat_angle(iC,azel_idx,1:2);  % azimuth and elevation for
        call_all_azel{iB}(iC,:,:) = A.proc.mic_to_bat_angle(iC,:,1:2);  % azimuth and elevation for
    end
end

mic_loc = A.mic_loc;

% Mic numbers for each sides, to set valid range for call_loc
left_idx = [1,16,6,3,2,8];
top_idx = [13,15,9,11,14,12,10];
bottom_idx = [27,4,26,5,17,24,7];
right_idx = [18,21,23,22,19,20];
front_idx = [29,30,28,34,33,25,32,31];
y(1) = min(mic_loc(right_idx,2));
y(2) = max(mic_loc(left_idx,2));
z(1) = max(mic_loc(bottom_idx,3));
z(2) = min(mic_loc(top_idx,3));
x(1) = 0.15;
x(2) = min(mic_loc(front_idx,1));

% Select only call_loc within range
call_loc = cell2mat(call_loc);  % collapse all call_loc together
call_max_azel = cell2mat(call_max_azel);  % beam aim in bat orientation
call_all_azel = cell2mat(call_all_azel);  % mic angle in bat orientation
in_range_idx = call_loc(:,1)>x(1) & call_loc(:,1)<x(2) &...  % select those within ranage
               call_loc(:,2)>y(1) & call_loc(:,2)<y(2) &...
               call_loc(:,3)>z(1) & call_loc(:,3)<z(2);
call_loc_sel = call_loc(in_range_idx,:);
call_max_azel_sel = call_max_azel(in_range_idx,:);
call_all_azel_sel = call_all_azel(in_range_idx,:,:);


% Plot selected call location
fig_track = figure;
corder = get(gca,'colororder');
hmic = plot3(mic_loc(:,1),mic_loc(:,2),mic_loc(:,3),'ko');
grid on
axis equal
hold on
for iB=1:size(track,1)
    tr = squeeze(track(iB,:,:));
    htrack = plot3(tr(:,1),tr(:,2),tr(:,3),'color',corder(1,:));
end
hcall = plot3(call_loc(:,1),call_loc(:,2),call_loc(:,3),'.','color',corder(2,:));
hcall_sel = plot3(call_loc_sel(:,1),call_loc_sel(:,2),call_loc_sel(:,3),'o','color',corder(3,:));
legend([hmic,htrack,hcall,hcall_sel],{'Mic','Bat track','All call location','Selected call location'},...
       'location','eastoutside');
saveas(gca,fullfile(save_path,sprintf('bat_track_call_loc_.fig')));
saveas(gca,fullfile(save_path,sprintf('bat_track_call_loc_.png')));

% Set bp param
bp_info.c = 344;  % sound speed [m/s]
bp_info.freq= 35e3;  % [Hz]
bp_type = 'piston';
if strcmp(bp_type,'piston');
    bp_info.type = 'piston';
    bp_info.a = 4e-3;  % aperture diameter [m]
    bp_info.k = 2*pi*bp_info.freq/bp_info.c;  % wavenumber
elseif strcmp(bp_type,'gaussian')
    bp_info.type = 'gaussian';
    bp_info.mu = 0;  % mean
    bp_info.sigma = 1;  % variance
end

% Model mic output
mic_dB = nan(size(call_all_azel_sel,1),size(call_all_azel_sel,2));
for iC = 1:size(call_max_azel_sel,1)
    mic_dB(iC,:) = model_beam(bp_info,call_max_azel_sel(iC,:),squeeze(call_all_azel_sel(iC,:,:)));
end

% Interpolation for beampattern surface
fig = figure;

for iC = 1:length(mic_dB)
rr = mic_dB(iC,:);
az = squeeze(call_all_azel_sel(iC,:,1))';
el = squeeze(call_all_azel_sel(iC,:,2))';

good_call = isgoodcall(az/pi*180,el/pi*180,rr,goodcall_angle_range,mic_area_frac);  % determine if good call

if good_call
    fprintf('Call#%02d is a good call\n',iC);
    maxref = max(rr);
    [azq,elq] = meshgrid(min(az):pi/180:max(az),min(el):pi/180:max(el));
    rrq = griddata(az,el,rr,azq,elq,'natural');
    rrq_norm = rrq-maxref;
    
%     vq_rbf = rbfinterp([azq(:)';elq(:)'],rbfcreate([az(:)';el(:)'],rr(:)','RBFFunction','multiquadrics'));
    nidx = ~isnan(rr);
    vq_rbf = rbfinterp([azq(:)';elq(:)'],rbfcreate([az(nidx)';el(nidx)'],rr(nidx),'RBFFunction','multiquadrics'));
    vq_rbf = reshape(vq_rbf,size(azq));
    vq_rbf_norm = vq_rbf-maxref;
    
    contour_vec_space = ceil(min([rrq_norm(~isinf(rrq_norm)&~isnan(rrq_norm));vq_rbf_norm(~isinf(vq_rbf_norm)&~isnan(vq_rbf_norm))])/3)*3:3:0;
    
    % Determine boundary of sampled area
    kk = boundary(az,el,0);  % outer boundary of all measured points
    in = inpolygon(azq,elq,az(kk),el(kk));
    
    % Find -3dB contour
    c3db_rrq = contourc(azq(1,:)/pi*180,elq(:,1)/pi*180,rrq_norm,[0 -3]);
    c_xy_rrq = c3db_rrq(:,2:c3db_rrq(2,1)-1)';
    c3db_vq = contourc(azq(1,:)/pi*180,elq(:,1)/pi*180,vq_rbf_norm,[0 -3]);
    c_xy_vq = c3db_vq(:,2:c3db_vq(2,1)-1)';
    
    % Plot
    figure(fig)
    cla
    subplot(121)
    himg = imagesc(azq(1,:)/pi*180,elq(:,1)/pi*180,rrq_norm);
    set(himg,'AlphaData',in);  % make NaN part transparent
    hold on
    contour(azq/pi*180,elq/pi*180,rrq_norm,contour_vec_space,'w');
    plot(az/pi*180,el/pi*180,'wx','linewidth',2);
    plot(c_xy_rrq(:,1),c_xy_rrq(:,2),'w','linewidth',2);
    caxis([-15 0])
    colorbar('location','southoutside')
    axis equal
    grid on
    hold off
    axis([-180 180 -90 90]);
    title('Natural neighbor');
        
    subplot(122)
    himg = imagesc(azq(1,:)/pi*180,elq(:,1)/pi*180,vq_rbf_norm);
    set(himg,'AlphaData',in);  % make NaN part transparent
    hold on
    contour(azq/pi*180,elq/pi*180,vq_rbf_norm,contour_vec_space,'w');
    plot(az/pi*180,el/pi*180,'wx','linewidth',2);
    plot(c_xy_vq(:,1),c_xy_vq(:,2),'w','linewidth',2);
    caxis([-15 0])
    colorbar('location','southoutside')
    axis equal
    grid on
    axis([-180 180 -90 90]);
    title('Radial basis function');
    
    suptitle(sprintf('Call #%02d',iC));
    % pause(0.5)
   
    % Fit ellipse
    try
        [zg, ag, bg, alphag] = fitellipse(c_xy_vq);
        plotellipse(zg, ag, bg, alphag, 'k');
        AR(iC) = ag/bg;
    catch
        disp(['Call #',num2str(iC),' cannot be fitted with ellipse']);
    end
    
    hold off
    saveas(fig,fullfile(save_path,sprintf('call_%02d.png',iC)));

else
    fprintf('Call#%02d is not a good call\n',iC);
    AR(iC) = NaN;
end

end

figure;
[count,bin] = hist(AR(~isnan(AR)&AR<2&AR>0),0:0.05:2);
hh = bar(bin,count/sum(count));
set(hh,'facecolor',corder(1,:))
xlim([bin(1) bin(end)])
ylim([0 0.35])
xlabel('Aspect ratio')
ylabel('Count')
grid
if strcmp(bp_type,'piston')
    title(sprintf('%s, radius %0.2f mm, freq %0.2f kHz',regexprep(bp_info.type,'^.','${upper($0)}'),bp_info.a*1e3,bp_info.freq/1e3));
    saveas(gca,fullfile(save_path,sprintf('%s_a%0.2fmm_%0.2fkHz_ARhist.fig',regexprep(bp_info.type,'^.','${upper($0)}'),bp_info.a*1e3,bp_info.freq/1e3)));
    saveas(gca,fullfile(save_path,sprintf('%s_a%0.2fmm_%0.2fkHz_ARhist.png',regexprep(bp_info.type,'^.','${upper($0)}'),bp_info.a*1e3,bp_info.freq/1e3)));
elseif strcmp(bp_type,'gaussian')
    title(sprintf('%s, mu %0.2f rad, sigma %0.2f rad',regexprep(bp_info.type,'^.','${upper($0)}'),bp_info.mu,bp_info.sigma));
    saveas(gca,fullfile(save_path,sprintf('%s_mu_%0.2frad_sigma%0.2frad.fig',regexprep(bp_info.type,'^.','${upper($0)}'),bp_info.mu,bp_info.sigma)));
    saveas(gca,fullfile(save_path,sprintf('%s_mu_%0.2frad_sigma%0.2frad.png',regexprep(bp_info.type,'^.','${upper($0)}'),bp_info.mu,bp_info.sigma)));
end
