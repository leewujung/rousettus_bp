% 2015 12 04  Plot all bat tracks and good call location

clear
usrn = getenv('username');

% Bat data path
base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
bat_proc_dir = './proc_output_rousettus_checked';
bat_proc_file = dir(fullfile(base_path,bat_proc_dir,'rousettus_20150825_*.mat'));
freq = 35e3;

iB = 10;
iC = 17;
A = load(fullfile(base_path,bat_proc_dir,bat_proc_file(iB).name));
track = A.track.track_interp;
% good_call_idx = find(A.proc.chk_good_call==1);
bat_loc = A.proc.bat_loc_at_call(iC,:);
call_loc_idx_on_track = A.proc.call_loc_on_track_interp(iC);
angle = squeeze(A.proc.mic_to_bat_angle(iC,:,1:2))/pi*180;
angle(:,1) = -angle(:,1);
mic_loc = A.mic_loc;
aim_v = A.head_aim.head_aim_int(call_loc_idx_on_track,:);  % head aim at call
norm_v = A.head_normal.head_normal_int(call_loc_idx_on_track,:);  % head normal at call

% Mic numbers for each sides
left_idx = [1,16,6,3,2,8];
top_idx = [13,15,9,11,14,12,10];
bottom_idx = [27,4,26,5,17,24,7];
right_idx = [18,21,23,22,19,20];
front_idx = [29,30,28,34,33,25,32,31];

% Spherical surface
r_proj = 0.7;
phi_lim = sort(angle(:,1));
phi_lim = [floor(phi_lim(1)/10)*10,ceil(phi_lim(end)/10)*10];
theta_lim = sort(angle(:,2));
theta_lim = [floor(theta_lim(1)/10)*10,ceil(theta_lim(end)/10)*10];

[phi,theta] = meshgrid(phi_lim(1):10:phi_lim(end),theta_lim(1):10:theta_lim(end));
[grid_x,grid_y,grid_z] = sph2cart(phi/180*pi,theta/180*pi,r_proj);
[mic_x,mic_y,mic_z] = sph2cart(angle(:,1)/180*pi,angle(:,2)/180*pi,r_proj);



R_prox = [1,0,0;...  % first rotation, from global to bat frame
          0,1,0;...
          0,0,1]';
xaxis_v = cross(aim_v,norm_v);  % x-axis (head aim x head normal)
xaxis_v = xaxis_v/norm(xaxis_v);  % normalized to unit vec
zaxis_v = cross(xaxis_v,aim_v);  % z-axis (x-axis x head aim)
R_dist = [xaxis_v;...
          aim_v;...
          zaxis_v]';
T1 = R_dist'*R_prox;

R_dist = [0,1,0;...  % second rotation, from bat aiming at +Y to bat aiming at +X
          -1,0,0;...
          0,0,1]';
T2 = R_dist'*R_prox;

TT = T2*T1;
TT_inv = inv(TT);

grid_inv = (TT_inv*([grid_x(:),grid_y(:),grid_z(:)]'))';
x_grid_inv = reshape(grid_inv(:,1),size(grid_x));
y_grid_inv = reshape(grid_inv(:,2),size(grid_y));
z_grid_inv = reshape(grid_inv(:,3),size(grid_z));
mic_inv = (TT_inv*([mic_x(:),mic_y(:),mic_z(:)]'))';
x_mic_inv = reshape(mic_inv(:,1),size(mic_x));
y_mic_inv = reshape(mic_inv(:,2),size(mic_y));
z_mic_inv = reshape(mic_inv(:,3),size(mic_z));

% Plot
mic_color = [1 1 1]*200/255;
sph_color = [0 1 1]*220/255;
bat_color = [0 0 0];

figure;
corder = get(gca,'colororder');
hmic_top = plot3(mic_loc(top_idx,1),mic_loc(top_idx,2),mic_loc(top_idx,3),'o','markerfacecolor',mic_color,'markersize',4,'markeredgecolor','none');
hold on
hmic_bottom = plot3(mic_loc(bottom_idx,1),mic_loc(bottom_idx,2),mic_loc(bottom_idx,3),'o','markerfacecolor',mic_color,'markersize',4,'markeredgecolor','none');
hmic_left = plot3(mic_loc(left_idx,1),mic_loc(left_idx,2),mic_loc(left_idx,3),'o','markerfacecolor',mic_color,'markersize',4,'markeredgecolor','none');
hmic_right = plot3(mic_loc(right_idx,1),mic_loc(right_idx,2),mic_loc(right_idx,3),'o','markerfacecolor',mic_color,'markersize',4,'markeredgecolor','none');
hmic_front = plot3(mic_loc(front_idx,1),mic_loc(front_idx,2),mic_loc(front_idx,3),'o','markerfacecolor',mic_color,'markersize',4,'markeredgecolor','none');
hs = surf(bat_loc(1)+x_grid_inv,bat_loc(2)+y_grid_inv,bat_loc(3)+z_grid_inv,'facecolor','none','edgecolor',sph_color);
hbat = plot3(bat_loc(1),bat_loc(2),bat_loc(3),'*','color',bat_color);
h_bataim = quiver3(bat_loc(1),bat_loc(2),bat_loc(3),aim_v(1),aim_v(2),aim_v(3),'k');
h_batnorm = quiver3(bat_loc(1),bat_loc(2),bat_loc(3),norm_v(1),norm_v(2),norm_v(3),'k');
h_batxdir = quiver3(bat_loc(1),bat_loc(2),bat_loc(3),xaxis_v(1),xaxis_v(2),xaxis_v(3),'k');
hmic_sph = plot3(bat_loc(1)+x_mic_inv,bat_loc(2)+y_mic_inv,bat_loc(3)+z_mic_inv,'o','markerfacecolor','r','markersize',4,'markeredgecolor','k');
axis equal
grid on
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');