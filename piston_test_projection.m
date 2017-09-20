% 2015 11 02  Use piston model to test projection onto az-el plane
% 2015 11 04  Make beam axis independent of bat head aim

clear
usrn = getenv('username');
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\beampattern_processing']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\fitellipse']);

% Param
mu = 0;
sigma = 0.5;
bp = 'piston';  % gauss, piston
aim_v = [0,1,0];
norm_v = [0,0,1];
beam_v = [1,1,0];
c = 344;
freq = 35e3;
a = 5e-3;

% Param proc
k = 2*pi*freq/c;
aim_v = aim_v/norm(aim_v);  % conver to unit vector
norm_v = norm_v/norm(norm_v);
beam_v = beam_v/norm(beam_v);

% Finely sampled beampattern
tt_fine = (-pi/2:pi/50:pi/2)';
pp_fine = (0:pi/50:pi)';
[t_fine,p_fine] = meshgrid(tt_fine,pp_fine);
t_fine = reshape(t_fine,1,[])';
p_fine = reshape(p_fine,1,[])';
[XXfine,YYfine,ZZfine] = sph2cart(p_fine,t_fine,ones(numel(p_fine),1));
pol_angle_fine = acos([XXfine,YYfine,ZZfine]*beam_v');
switch bp
    case 'piston'
        rr_fine = 20*log10(abs(2*besselj(1,k*a*sin(pol_angle_fine))./(k*a*sin(pol_angle_fine))));
        rr_fine(pol_angle_fine==0) = 0;
    case 'gauss'
        rr_fine = 20*log10(normpdf(pol_angle_fine,mu,sigma));
end
figure;
mesh(pp_fine,tt_fine,reshape(rr_fine,length(pp_fine),[])');


% Randomly sampled beampattern
theta = (linspace(0,1,6)-0.5)'*pi;
phi = linspace(0,1,6)'*pi;

[theta,phi] = meshgrid(theta,phi);
theta = reshape(theta,1,[])';
phi = reshape(phi,1,[])';

theta = theta+(rand(length(theta),1)-0.5)*0.2*pi;
phi = phi+(rand(length(phi),1)-0.5)*0.2*pi;

% theta = (rand(32,1)-0.5)*pi;
% phi = rand(32,1)*pi;

% theta = (rand(100,1)-0.5)*pi;
% phi = rand(100,1)*pi;


% Calculate beampattern
[XX,YY,ZZ] = sph2cart(phi,theta,ones(numel(phi),1));
pol_angle = acos([XX,YY,ZZ]*beam_v');
switch bp
    case 'piston'
        rr = 20*log10(abs(2*besselj(1,k*a*sin(pol_angle))./(k*a*sin(pol_angle))));
        rr(pol_angle==0) = 0;
    case 'gauss'
        rr = 20*log10(normpdf(pol_angle,mu,sigma));
end

% Calculate angle to bat
[mic2bat_2d,~] = find_mic_az_el_to_bat_fcn([XX,YY,ZZ],aim_v,norm_v);
az = mic2bat_2d(:,1);
el = mic2bat_2d(:,2);

ref = max(max(rr));
[azq,elq] = meshgrid(min(az):pi/180:max(az),min(el):pi/180:max(el));
rrq = griddata(az,el,rr,azq,elq,'natural');
rrq_norm = rrq-ref;

vq_rbf = rbfinterp([azq(:)';elq(:)'],rbfcreate([az(:)';el(:)'],rr(:)','RBFFunction','multiquadrics'));
vq_rbf = reshape(vq_rbf,size(azq));
vq_rbf_norm = vq_rbf-ref;

contour_vec_space = ceil(min([rrq_norm(~isinf(rrq_norm)&~isnan(rrq_norm));vq_rbf_norm(~isinf(vq_rbf_norm)&~isnan(vq_rbf_norm))])/3)*3:3:0;

% Determine boundary of sampled area
kk = boundary(az,el,0);  % outer boundary of all measured points
in = inpolygon(azq,elq,az(kk),el(kk));
c3db_rrq = contourc(azq(1,:)/pi*180,elq(:,1)/pi*180,rrq_norm,[0 -3]);
c_xy_rrq = c3db_rrq(:,2:c3db_rrq(2,1)-1)';
c3db_vq = contourc(azq(1,:)/pi*180,elq(:,1)/pi*180,vq_rbf_norm,[0 -3]);
c_xy_vq = c3db_vq(:,2:c3db_vq(2,1)-1)';

% Plot
figure;
subplot(121)
himg = imagesc(azq(1,:)/pi*180,elq(:,1)/pi*180,rrq_norm);
set(himg,'AlphaData',in);  % make NaN part transparent
hold on
contour(azq/pi*180,elq/pi*180,rrq_norm,contour_vec_space,'w');
plot(az/pi*180,el/pi*180,'wx','linewidth',2);
plot(c_xy_rrq(:,1),c_xy_rrq(:,2),'w','linewidth',2);
caxis([-30 0])
colorbar('location','southoutside')
axis equal
grid on

subplot(122)
himg = imagesc(azq(1,:)/pi*180,elq(:,1)/pi*180,vq_rbf_norm);
set(himg,'AlphaData',in);  % make NaN part transparent
hold on
contour(azq/pi*180,elq/pi*180,vq_rbf_norm,contour_vec_space,'w');
plot(az/pi*180,el/pi*180,'wx','linewidth',2);
plot(c_xy_vq(:,1),c_xy_vq(:,2),'w','linewidth',2);
caxis([-30 0])
colorbar('location','southoutside')
axis equal
grid on
% Fit ellipse
% sss = input('1-fit ellipse, 0-no fitting:');
% if sss
    [zg, ag, bg, alphag] = fitellipse(c_xy_vq);
    plotellipse(zg, ag, bg, alphag, 'k')
% end

