% 2015 11 04  Test rotating max amplitude point to middle of az-el plane
% 2015 11 20  Use new format in call_dB selection

% clear
usrn = getenv('username');
DATA_DIR = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\proc_output'];
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);

bpfile = 'rousettus_20150825_36134_02_mic_data_bp_proc.mat';

data = load(fullfile(DATA_DIR,bpfile));

plot_opt = 0;
freq_wanted = 35e3;
% angle_range = [-60 60];

iC = 17;

% Call info
call_dB = nan(data.mic_data.num_ch_in_file,1);
for iM=1:data.mic_data.num_ch_in_file
    freq = data.proc.call_freq_vec{iC,iM};
    [~,fidx] = min(abs(freq-freq_wanted));
    call_dB(iM) = data.proc.call_psd_dB_comp_re20uPa_withbp{iC,iM}(fidx);
end
% [~,fidx] = min(abs(freq-data.proc.call_freq_vec{iC}));
% call_dB = squeeze(data.proc.call_psd_dB_comp_re20uPa_withbp{iC}(fidx,:));
ch_good_loc = ~isnan(data.mic_loc(:,1))';
ch_ex_sig = find(isnan(call_dB));  % low quality channel from call extraction function
mic_num = 1:data.mic_data.num_ch_in_file;
angle_notnanidx = ~ismember(1:data.mic_data.num_ch_in_file,ch_ex_sig) & ch_good_loc;
az = squeeze(data.proc.mic_to_bat_angle(iC,angle_notnanidx,1))'/pi*180;
el = squeeze(data.proc.mic_to_bat_angle(iC,angle_notnanidx,2))'/pi*180;
bnd_idx = boundary(az,el,0);  % outer boundary of all measured points
bnd = [az(bnd_idx),el(bnd_idx)];

% Interpolation
[maxref,maxref_idx] = max(call_dB(angle_notnanidx));
[azq,elq] = meshgrid(min(az):max(az),min(el):max(el));
vq_n = griddata(az,el,call_dB(angle_notnanidx),azq,elq,'natural');  % natural neighbor interp
vq_n_norm = vq_n-maxref;
vq_r = rbfinterp([azq(:)';elq(:)'],rbfcreate([az(:)';el(:)'],call_dB(angle_notnanidx)','RBFFunction','multiquadrics'));
vq_r = reshape(vq_r,size(azq));
vq_r_norm = vq_r-maxref;
[in,on] = inpolygon(azq(:),elq(:),bnd(:,1),bnd(:,2));
in_smpl_poly = in|on;
clear in on
vq_r(~in_smpl_poly) = NaN;
vq_r_norm(~in_smpl_poly) = NaN;

% Rotation


% Plot
figure
axesm eckert4;
framem; gridm;
axis off
elqm = elq;
elqm(isnan(vq_r_norm)) = NaN;
azqm = azq;
azqm(isnan(vq_r_norm)) = NaN;
gg = geoshow(elqm,azqm,vq_r_norm,'displaytype','texturemap');



