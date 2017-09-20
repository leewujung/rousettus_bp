% 2015 10 25  Test RBF interpolation

username = getenv('username');
addpath(['C:\Users\',username,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
addpath(['C:\Users\',username,'\Dropbox\0_CODE\MATLAB\fitellipse']);

% Load A
A=load('C:\Users\Wu-Jung Lee\Dropbox\0_ANALYSIS\bp_processing\proc_output\rousettus_20150825_36134_02_mic_data_bp_proc.mat');

for iC=10:length(A.mic_data.call_idx_w_track)
% call_num = 17;
call_num = iC;
freq_wanted = 35e3;

mic_to_bat_angle = squeeze(A.proc.mic_to_bat_angle(call_num,:,:));
[~,fidx] = min(abs(freq_wanted-A.proc.call_freq_vec{call_num}));
call_dB = squeeze(A.proc.call_psd_dB_comp_re20uPa_withbp{call_num}(fidx,:));

% Check for channels to be excluded
if isempty(A.proc.ch_ex{call_num})
    ch_ex_manual = [];
else
    ch_ex_manual = A.proc.ch_ex{call_num};
end
ch_ex_sig = find(isnan(call_dB));  % low quality channel from call extraction function

% Check for mics without location A
ch_good_loc = ~isnan(A.mic_loc(:,1))';

% Matlab interpolation
contour_vec = 0:-3:-60;

mic_num = 1:A.mic_data.num_ch_in_file;
angle_notnanidx = ~ismember(1:A.mic_data.num_ch_in_file,union(ch_ex_manual,ch_ex_sig)) & ch_good_loc;
az = mic_to_bat_angle(angle_notnanidx,1);
el = mic_to_bat_angle(angle_notnanidx,2);

[azq,elq] = meshgrid(min(az):pi/180:max(az),min(el):pi/180:max(el));
v = call_dB(angle_notnanidx);
vq_linear = griddata(az,el,v,azq,elq,'linear');
vq_natural = griddata(az,el,v,azq,elq,'natural');
vq_cubic = griddata(az,el,v,azq,elq,'cubic');


% RBF interpolation
vq_rbf_gauss = rbfinterp([azq(:)';elq(:)'],rbfcreate([az(:)';el(:)'],v(:)','RBFFunction','gaussian','RBFConstant',0.5));
vq_rbf_gauss = reshape(vq_rbf_gauss,size(azq));
vq_rbf_cubic = rbfinterp([azq(:)';elq(:)'],rbfcreate([az(:)';el(:)'],v(:)','RBFFunction','cubic'));
vq_rbf_cubic = reshape(vq_rbf_cubic,size(azq));
vq_rbf_multiq = rbfinterp([azq(:)';elq(:)'],rbfcreate([az(:)';el(:)'],v(:)','RBFFunction','multiquadrics'));
vq_rbf_multiq = reshape(vq_rbf_multiq,size(azq));
option_multiq = rbfcreate([az(:)';el(:)'],v(:)','RBFFunction','multiquadrics');

% Determine polygon
k = boundary(az,el,0);  % outer boundary of all measured points
in = inpolygon(azq,elq,az(k),el(k));

% Test fitting ellipse
figure;
himg = imagesc(azq(1,:)/pi*180,elq(:,1)/pi*180,vq_rbf_multiq);
set(himg,'AlphaData',in);  % make NaN part transparent
hold on
[c,hc] = contour(azq/pi*180,elq/pi*180,vq_rbf_multiq-max(max(vq_rbf_multiq)),[0 -3],'w','linewidth',2);
c_xy = c(:,2:c(2,1)-1)';
plot(c_xy(:,1),c_xy(:,2),'w','linewidth',2);
axis xy
axis equal
title(['Call #',num2str(iC)]);

% Fit elliipse
sss = input('1-fit ellipse, 0-no fitting:');
if sss
    [zg, ag, bg, alphag] = fitellipse(c_xy);
    plotellipse(zg, ag, bg, alphag, 'k')
    pause
end

% 
% % Comparison
% figure('units','normalized','outerposition',[0 0 1 1])
% for iF=1:6
%     switch iF
%         case 1
%             vvv = vq_linear;
%             ttt = 'Matlab griddata, linear';
%         case 3
%             vvv = vq_natural;
%             ttt = 'Matlab griddata, natural neighbor';
%         case 5
%             vvv = vq_cubic;
%             ttt = 'Matlab griddata, cubic';
%         case 2
%             vvv = vq_rbf_gauss;
%             ttt = 'RBF, gaussian';
%         case 4
%             vvv = vq_rbf_multiq;
%             ttt = 'RBF, multiquadric';
%         case 6
%             vvv = vq_rbf_cubic;
%             ttt = 'RBF, cubic';
%     end
%     subplot(3,2,iF);
%     himg = imagesc(azq(1,:)/pi*180,elq(:,1)/pi*180,vvv);
%     set(himg,'AlphaData',in);  % make NaN part transparent
%     hold on
%     contour(azq/pi*180,elq/pi*180,vvv-max(max(vvv)),[0 -3],'w','linewidth',2);
%     plot(az/pi*180,el/pi*180,'kx');
%     title(ttt);
%     axis xy
%     axis equal
% end
% suptitle(['Call #',num2str(iC)]);

% pause
end
