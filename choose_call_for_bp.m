% 2015 11 04  Choose which calls to be included in bp analysis

% clear
usrn = getenv('username');
DATA_DIR = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\proc_output'];
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);

bpfile = 'rousettus_20150825_36134_03_mic_data_bp_proc.mat';

A = load(fullfile(DATA_DIR,bpfile));

plot_opt = 0;
freq = 35e3;
angle_range = [-60 60];

for iC = 18%1:length(A.mic_data.call_idx_w_track)
    
% Get valid mic channels and proc output
[~,fidx] = min(abs(freq-A.proc.call_freq_vec{iC}));
call_dB = squeeze(A.proc.call_psd_dB_comp_re20uPa_withbp{iC}(fidx,:));
ch_good_loc = ~isnan(A.mic_loc(:,1))';
ch_ex_sig = find(isnan(call_dB));  % low quality channel from call extraction function
mic_num = 1:A.mic_data.num_ch_in_file;
angle_notnanidx = ~ismember(1:A.mic_data.num_ch_in_file,ch_ex_sig) & ch_good_loc;
az = squeeze(A.proc.mic_to_bat_angle(iC,angle_notnanidx,1))'/pi*180;
el = squeeze(A.proc.mic_to_bat_angle(iC,angle_notnanidx,2))'/pi*180;
bnd_idx = boundary(az,el,0);  % outer boundary of all measured points
bnd = [az(bnd_idx),el(bnd_idx)];

% Interpolation
[maxref,maxref_idx] = max(call_dB(angle_notnanidx));
[azq,elq] = meshgrid(min(az):max(az),min(el):max(el));
vq_n = griddata(az,el,call_dB(angle_notnanidx),azq,elq,'natural');  % natural neighbor interp
vq_n_norm = vq_n-maxref;
vq_r = rbfinterp([azq(:)';elq(:)'],rbfcreate([az(:)';el(:)'],call_dB(angle_notnanidx),'RBFFunction','multiquadrics'));
vq_r = reshape(vq_r,size(azq));
vq_r_norm = vq_r-maxref;
[in,on] = inpolygon(azq(:),elq(:),bnd(:,1),bnd(:,2));
in_smpl_poly = in|on;
clear in on
vq_r(~in_smpl_poly) = NaN;
vq_r_norm(~in_smpl_poly) = NaN;

% % 3dB contour
% c3db_vq_n = contourc(azq(1,:),elq(:,1),vq_n_norm,[0 -3]);
% c3db_vq_n = c3db_vq_n(:,2:c3db_vq_n(2,1)-1)';
% c3db_vq_r = contourc(azq(1,:),elq(:,1),vq_r_norm,[0 -3]);
% c3db_vq_r = c3db_vq_r(:,2:c3db_vq_r(2,1)-1)';
% c3db_vq_r_in_3db_idx = inpolygon(c3db_vq_r(:,1),c3db_vq_r(:,2),bnd(:,1),bnd(:,2));
% c3db_vq_r = c3db_vq_r(c3db_vq_r_in_3db_idx,:);  % only take the 3dB contour within the sample range

% % more than 10 mic within az/el range
% azel_in_range = az<angle_range(2) & az>angle_range(1) &...
%                 el<angle_range(2) & el>angle_range(1);
% angle_flag = sum(azel_in_range)>=10;

% more than 10 mic within az/el range around peak mic
azel_in_range = az<az(maxref_idx)+angle_range(2) & az>az(maxref_idx)+angle_range(1) &...
                el<el(maxref_idx)+angle_range(2) & el>el(maxref_idx)+angle_range(1);
angle_flag = sum(azel_in_range)>=10;

% highest amp within az/el range
[~,max_amp_idx] = max(A.proc.call_psd_dB_comp_re20uPa_withbp{iC}(fidx,:));
max_angle_flag = az(max_amp_idx)<angle_range(2) & az(max_amp_idx)>angle_range(1) &...
                 el(max_amp_idx)<angle_range(2) & el(max_amp_idx)>angle_range(1);

% highest amp not on edge
max_edge_flag = ~ismember(max_amp_idx,bnd_idx);

% at least 3 mics within 3dB contour
larger_3db_idx = (call_dB-maxref)>-3.5;
ch_num_3db_flag = sum(larger_3db_idx)>=3;

% Good call?
good_call(iC) = angle_flag & max_angle_flag & max_edge_flag & ch_num_3db_flag;

% Plot
if plot_opt==1
    azq_unique = unique(azq);
    elq_unique = unique(elq);
    scrsz = get(groot,'ScreenSize');
    figure('Position',[500 150 800 700]);
    subplot(211)
    himg = contourf(azq,elq,vq_r_norm,-3:-3:-180,'w');
    hold on
    plot(az,el,'mx','linewidth',1);
    contour(azq,elq,vq_r_norm,-3:-3:-180,'w');
    axis xy
    axis equal
    axis([-180 180 -90 90]);
    title(sprintf('Call#%d on track',iC));
    if good_call==1
        text(-170,-70,'GOOD','fontsize',18,'color','m');
    else
        text(-170,-70,'BAD','fontsize',18,'color','m');
    end
    subplot(212)
    himg = imagesc(azq_unique,elq_unique,vq_r_norm);
    set(himg,'alphadata',~isnan(vq_n));
    hold on
    plot(az,el,'mx','linewidth',1);
    contour(azq,elq,vq_r_norm,-3:-3:-180,'w');
    axis xy
    axis equal
    axis([-180 180 -90 90]);
    title(sprintf('Call#%d on track',iC));
    if good_call==1
        text(-170,-70,'GOOD','fontsize',18,'color','m');
    else
        text(-170,-70,'BAD','fontsize',18,'color','m');
    end
end
% pause
% close
end

