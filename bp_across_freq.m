% 2015 11 28  3D beampattern across differen freq

clear
warning off

usrn = getenv('username');
base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
bat_proc_path = './proc_output';
bat_proc_file = 'rousettus_20150825_36134_02_mic_data_bp_proc';

addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\beampattern_processing']);

data = load(fullfile(base_path,bat_proc_path,bat_proc_file));
good_call_idx = find(data.proc.chk_good_call);
call_seq = 22;
 
% Get angle & good mic_idx info
[~,az,el,ch_include_idx] = get_call_azel_dB_data(data,35e3,call_seq);

freq_wanted = (20:5:50)*1e3;
for iM=1:data.mic_data.num_ch_in_file
    % Interpolate across freq
    dB = data.proc.call_psd_dB_comp_re20uPa_withbp{call_seq,iM};
    freq_meas = data.proc.call_freq_vec{call_seq,iM};
    call_dB(:,iM) = interp1(freq_meas,dB,freq_wanted);
end

num_freq = length(freq_wanted);
for iF=1:num_freq
    fprintf('Frequency: %2.1f kHz\n',freq_wanted(iF)/1e3);
    
    % Interpolate across mics
    [vq,vq_norm,azq,elq] = interp_bp(az(ch_include_idx),el(ch_include_idx),call_dB(iF,ch_include_idx),'rbf');
    
    % Project lat-lon to map projection distance
    mstruct = defaultm('eckert4');
    mstruct = defaultm(mstruct);
    [xq,yq] = mfwdtran(mstruct,elq/pi*180,azq/pi*180);

    % Find 3dB contour
    figure
    [C,h] = contour(xq,yq,vq_norm,0:-3:-45,'fill','on');
    Cout = parse_contour_output(C);
    c3db_xy = [];
    for iT=1:length(Cout)  % in case contour break into pieces
        if Cout(iT).Level == -3
            c3db_xy = [c3db_xy; Cout(iT).X',Cout(iT).Y'];
        end
    end
%     pause
    close
    c3db_all{iF} = c3db_xy;
end

[xqlim,yqlim] = mfwdtran(mstruct,[-90 90 0 0],[0 0 -180 180]);
xqlim = xqlim(3:4);
yqlim = yqlim(1:2);

colorset = jet(num_freq);
 
fig_3db = figure;
hold on
for iF=1:num_freq
    plot3(c3db_all{iF}(:,1),c3db_all{iF}(:,2),repmat(freq_wanted(iF)/1e3*0.01,size(c3db_all{iF},1)),...
          'linewidth',2,'color',colorset(iF,:));
end
axis equal
title(sprintf('Call# %02d',call_seq))
colormap(jet(num_freq))
colorbar('Ticks',linspace(0,1,num_freq),'TickLabels',{num2str(freq_wanted'/1e3)})
axis([xqlim,yqlim])
grid

% c3db_clpse = cell2mat(c3db_all');
% xx3 = min(c3db_clpse(:,1)):range(c3db_clpse(:,1))/100:max(c3db_clpse(:,1));
% yy3 = min(c3db_clpse(:,2)):range(c3db_clpse(:,2))/100:max(c3db_clpse(:,2));
% [xq3,yq3] = meshgrid(xx3,yy3);
% V = zeros(size(xq3,1),size(xq3,2),num_freq);
% for iF=1:num_freq
%     [in,on] = inpolygon(xq3,yq3,c3db_all{iF}(:,1),c3db_all{iF}(:,2));
%     in = in|on;
%     Vtmp = zeros(size(xq3));
%     Vtmp(in) = 1;
%     V(:,:,iF) = reshape(Vtmp,size(xq3));
% end


