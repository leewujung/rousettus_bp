% 2016 08 18  Plot piston model bp and multi-freq contours
% 2017 04 12  New piston model fitting using new results with corrected bat_to_mic_angle

clear
usrn = getenv('username');
% Set up various paths
base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp\analysis_results_figs';
addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
addpath('F:\Dropbox\0_CODE\OpenBEM_02-2015\3DBEM\testcases');


[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

save_fname = script_name;

%% Globe bp plot
freq = 35e3;
c = 344;
a = 4.91*1e-3; % [m]  % before 2017: a=4.92*1e-3
k = 2*pi*freq/c;

[theta,phi] = meshgrid(-90:5:90,-180:5:180);  % field point locations
sz = size(theta);
theta = theta(:);
phi = phi(:);
[x,y,z] = sph2cart(phi/180*pi,theta/180*pi,1);

vec0 = [1,0,0];
angle = acos([x,y,z]*vec0');

angle = reshape(angle,sz);
theta = reshape(theta,sz);
phi = reshape(phi,sz);

piston_bp = 20*log10(abs(2*besselj(1,k*a*sin(angle))./...
    (k*a*sin(angle))));
piston_bp(isnan(piston_bp)) = 0;
piston_bp(angle>pi/2) = min(piston_bp(:));

% plot
map_proj = 'eckert4';
az_plot_limit = 180;
cvec = 0:-3:-30;
vq_norm_min = -27;
contour_vec = 0:-3:(floor(vq_norm_min/3)-1)*3;
cvec_min_idx = find(contour_vec-vq_norm_min<0,1,'first');
        
M.el = theta;
M.az = phi;
M.pp_plot = piston_bp;

fig_bp = figure;
plot_model_bp_globe(gca,M,cvec,map_proj,az_plot_limit);  % bp globe
colorbar('southoutside','ticks',fliplr(contour_vec(1:cvec_min_idx)));
colormap(parula(cvec_min_idx-1));
caxis([contour_vec(cvec_min_idx) 0]);
tightmap
        
saveSameSize(fig_bp,'file',fullfile(save_path,[save_fname,'_bp.png']),...
    'format','png','renderer','painters');
epswrite(fullfile(save_path,[save_fname,'_bp.eps']));


%% Multi-freq contours
freq_all = (25:5:55)*1e3;  % [kHz]
angle = 0:0.0001*pi:pi/2;

mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);
cgrey = 200*ones(1,3)/255;
num_freq = length(freq_all);

fig_cntr = figure;
axesm(map_proj);
axis off
gridm('gcolor',cgrey,'glinestyle','-');
framem('fedgecolor',cgrey,'flonlimit',[-1 1]*az_plot_limit,...
    'flinewidth',1);
colorset = jet(num_freq);

for iF=1:num_freq
    k = 2*pi*freq_all(iF)/c;

    piston_bp = 20*log10(abs(2*besselj(1,k*a*sin(angle))./...
        (k*a*sin(angle))));
    [~,idx_3db] = min(abs(piston_bp+3));
    angle_sel = angle(idx_3db);
    
    theta = angle_sel*sin(0:0.001*pi:2*pi);
    phi = angle_sel*cos(0:0.001*pi:2*pi);
    
    % plot
    figure(fig_cntr)
    plotm(theta/pi*180,phi/pi*180,'color',colorset(iF,:),...
        'linewidth',2);
    hold on
    
end
colormap(jet(num_freq))
colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
    'TickLabels',{num2str(freq_all'/1e3)},'location','southoutside');
tightmap
set(gcf,'color','w')

saveSameSize(fig_cntr,'file',fullfile(save_path,[save_fname,'_cntr.png']),...
    'format','png','renderer','painters');
epswrite(fullfile(save_path,[save_fname,'_cntr.eps']));










