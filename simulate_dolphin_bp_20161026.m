% 2016 10 26  Make up dolphin path, click direction, and mic location

clear
usrn = getenv('username');
% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp\analysis_results_figs';
    addpath('F:\Dropbox\0_CODE\beampattern_other_code');
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
else
    base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp\analysis_results_figs'];
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\beampattern_other_code']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
end

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

mic_spacing = 0.25;
screen_size = [2,2];  % x/y dimension of screen [m]
rep_num = 4;

S.mic_spacing = mic_spacing;
S.screen_size = screen_size;

save_fname = sprintf('%s_micsp%2.2fm',script_name,mic_spacing);


% % mic_config = 1;  % 0-random, 1-offset rows
% mic_num = 16;

% Seed rng
rng_seed = 168;
rng(rng_seed);  % seed the random number generator

S.rng_seed = rng_seed;


% Mic locs
% mic_x_temp = linspace(0.15,1.75,floor(1.7/mic_spacing));  % 0.4 spacing
mic_x_temp = linspace(0.15,1.85,floor(1.7/mic_spacing));  % 0.25 spacing
mic_x = repmat(mic_x_temp,length(mic_x_temp),1);
idx_odd = 2:2:length(mic_x_temp);
for iO=1:length(idx_odd)
    mic_x(idx_odd(iO),:) = mic_x_temp+mic_spacing/2;
end
mic_y = repmat(mic_x_temp',1,length(mic_x_temp));
mic_sz = size(mic_x);
mic_loc = [mic_x(:),mic_y(:)];
mic_loc = [mic_loc,zeros(size(mic_loc,1),1)];
[xq,yq] = meshgrid(0:0.01:2,0:0.01:2);

S.mic_loc = mic_loc;


% Dolphin locs
d_num = 15;
dloc_bnd = [2,2,4];
d_loc(:,1) = rand(d_num,1)*dloc_bnd(1);
d_loc(:,2) = rand(d_num,1)*dloc_bnd(2);
d_loc(:,3) = linspace(0.5,dloc_bnd(3),d_num);
% d_loc(:,3) = rand(d_num,1)*dloc_bnd(3)+0.5;
d_loc = flipud(sortrows(d_loc,3));
d_loc = repmat(d_loc,1,rep_num);  % rep rep_num times
d_loc = reshape(d_loc',3,[])';

S.d_num = d_num;
S.d_loc = d_loc;


% Click aim
c_num = d_num*rep_num;
cloc_bnd = [2,2];
c_loc(:,1) = rand(c_num,1)*cloc_bnd(1);
c_loc(:,2) = rand(c_num,1)*cloc_bnd(2);
c_loc(:,3) = zeros(size(c_loc,1),1);

S.c_num = c_num;
S.c_loc = c_loc;


% Plot setup
fig_setup = figure;
plot3(mic_loc(:,1),mic_loc(:,2),mic_loc(:,3),'.')
hold on
plot3(d_loc(:,1),d_loc(:,2),d_loc(:,3),'-o')
plot3(c_loc(:,1),c_loc(:,2),c_loc(:,3),'x')
axis equal
grid
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');

saveas(fig_setup,fullfile(save_path,[save_fname,'_setup.fig']),'fig');
saveSameSize(fig_setup,'file',fullfile(save_path,[save_fname,'_setup.png']),...
    'format','png','renderer','painters');

% Beampattern params
match_deg = 5;  % [deg]
bp_info.c = 1500;  % sound speed [m/s]
bp_info.freq = 60e3;  % [Hz]
bp_info.k = 2*pi*bp_info.freq/bp_info.c;  % wavenumber

% Find aperture with best-fitting -3dB contour to mean azimuth of data
a_fine = (2:0.01:7)*1e-2; % [m]
theta_fine = 0:pi/1000:pi/2;
piston_bp = 20*log10(abs(2*besselj(1,bp_info.k*a_fine'*sin(theta_fine))./...
    (bp_info.k*a_fine'*sin(theta_fine))));
[C,~] = contour(theta_fine/pi*180,a_fine,piston_bp,[-3 -6],'fill','on');
c_level = parse_contour_output(C);
piston_3dB_c = [c_level(2).X; c_level(2).Y]';

[~,idx_mean] = min(abs(piston_3dB_c(:,1)-match_deg));
a_mean = piston_3dB_c(idx_mean,2);
close
bp_info.a = a_mean;

piston_bp_best = 20*log10(abs(2*besselj(1,bp_info.k*a_mean*sin(theta_fine))./...
    (bp_info.k*a_mean*sin(theta_fine))));

% Plot params
cvec = [0:-1:-30];

% Beampattern simulation
fig_bp = figure;
corder = get(gca,'colororder');
for iD=1:d_num*rep_num
    % Mic to bat vector
    mic2bat = mic_loc-repmat(d_loc(iD,:),size(mic_loc,1),1);
    mic2bat_len = sqrt(diag(mic2bat*mic2bat'));
    mic2bat = mic2bat./repmat(mic2bat_len,1,3);
    
    % Beam aim to bat vector
    aim_v = c_loc(iD,:)-d_loc(iD,:);
    aim_v = aim_v/norm(aim_v);
    
    pol_angle = acos(diag(mic2bat*repmat(aim_v',1,size(mic2bat,1))));

    mic_dB = 20*log10(abs(2*besselj(1,bp_info.k*bp_info.a*sin(pol_angle))./...
        (bp_info.k*bp_info.a*sin(pol_angle))));
    mic_dB(pol_angle==0) = 0;
    
    % Interpolation
    vq = griddata(mic_loc(:,1),mic_loc(:,2),mic_dB,xq,yq,'natural');
    
    [c_level,c_level_nan] = get_main_contour(vq,unique(xq),unique(yq),-3);
    
    % Plot
    figure(fig_bp)
    cla
    imagesc(unique(xq),unique(yq),vq);
    hold on
    plot(mic_loc(:,1),mic_loc(:,2),'o','linewidth',2,'color','w');
    plot(c_loc(iD,1),c_loc(iD,2),'*','linewidth',2,'color','r','markersize',10);
    if ~isempty(c_level)
        plot(c_level(:,1),c_level(:,2),'r')
    end
    axis xy
    xlabel('X (m)');
    ylabel('Y (m)');
%     colormap(parula(length(cvec)-1))
    colorbar('ticks',[-50:10:-10,-3,0]);
%     colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
%     'TickLabels',{num2str(freq_all'/1e3)},'location','southoutside');
%     caxis([0,-3,-10,-20:10:-60])
    caxis([-50 0])
    axis equal
    axis([0 2 0 2])
    set(gca,'xtick',0:0.5:2,'ytick',0:0.5:2);
    title(sprintf('Click %02d, Dist=%2.2f m',iD,d_loc(iD,3)));
    
    save_fname_curr = sprintf('%s_c%02d_d%2.2fm_bp',save_fname,iD,d_loc(iD,3));
    saveSameSize_100(fig_bp,'file',fullfile(save_path,[save_fname_curr,'.png']),...
    'format','png','renderer','painters');
end









