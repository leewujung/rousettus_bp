% 2015 11 25  Different aperture model

%% Rectangular aperture
a = 6e-3;  % height
b = 2e-2;  % width
freq = 35e3;
c = 344;
k = 2*pi*freq/c;
[tt,pp] = meshgrid(-90:0.5:90,-90:0.5:90);  % [deg]
theta = tt(:)/180*pi;  % [rad]
phi = pp(:)/180*pi;
alpha = k*a/2*sin(theta);
beta = k*b/2*sin(phi);
I_rect = 10*log10(sinc(alpha).^2.*sinc(beta).^2);
I_rect = reshape(I_rect,size(tt));

% Project lat-lon to map projection distance
mstruct = defaultm('eckert4');
mstruct = defaultm(mstruct);
[xq,yq] = mfwdtran(mstruct,tt,pp);

figure
% contourf(xq,yq,I_rect,[-30:3:0],'edgecolor','none');
surf(xq,yq,I_rect,'edgecolor','none');
caxis([-30 0])
colormap(parula(9))
view([0 0 1]);
axis equal


%% Circular aperture with rotation
a0 = 10e-3;
freq = 35e3;
c = 344;
k = 2*pi*freq/c;
[pp,tt] = meshgrid(0:0.5:180,-180:0.5:180);  % [deg]
I_circ = 20*log10(abs(2*besselj(1,k*a0*sin(tt/180*pi))./...
         (k*a0*sin(tt/180*pi))));
I_circ(isnan(I_circ)) = 0;
tt = tt-90;  % shift theta due to coord system definition
[tt_rot,pp_rot] = rotatem(tt,pp,[90 0],'forward','degrees');
mstruct = defaultm('ortho');
mstruct = defaultm(mstruct);
[xq_rot,yq_rot] = mfwdtran(mstruct,tt_rot,pp_rot);

figure
% contourf(xq_rot,yq_rot,I_circ,[-30:3:0],'edgecolor','none');
surf(xq_rot,yq_rot,I_circ,'edgecolor','none');
caxis([-30 0])
colormap(parula(9))
view([0 0 1]);
axis equal


%% Linear array with multiple elements (identical, isotropic elements)
d = 0.002;
n = 5;
freq = [20:5:50]*1e3;
c = 344;
k = 2*pi*freq/c;
% delta = d/4*k;
delta = 0.5;
% delta = 0;
theta = -90:90;
% psi = bsxfun(@times,k',d*sin(repmat(theta/180*pi,length(delta),1)+repmat(delta',1,length(theta))));
psi = bsxfun(@times,k',d*sin(repmat(theta/180*pi,length(delta),1)))+repmat(delta,length(k),length(theta));
E = sin(n*psi/2)./sin(psi/2)/n;
E(isnan(E)) = 1;
E = 20*log10(abs(E));
[~,mmidx] = max(E,[],2);
imagesc(theta,freq/1e3,E);
caxis([-30 0])
colormap(parula(10))
colorbar

%% Linear array with multiple elements (identical, isotropic elements)
d = 0.002;
n = 5;
freq = [20:2:50]*1e3;
c = 344;
k = 2*pi*freq/c;
delta = d/4*k;
% delta = 0;
theta = -90:90;
psi = bsxfun(@times,k',d*sin(repmat(theta/180*pi,length(delta),1)+repmat(delta',1,length(theta))));
E = sin(n*psi/2)./sin(psi/2)/n;
E(isnan(E)) = 1;
E = 20*log10(abs(E));
[~,mmidx] = max(E,[],2);


%% Linear array, direct summation with approximation
d = -0.002;
n = 5;
freq = 30*1e3;
c = 344;
theta = 0:0.5:180;
k = 2*pi*freq/c;
D = (0:n-1)*d;
b = sum(exp(-1j*k*D'*cos(theta/180*pi)),1);
B = 20*log10(abs(b));
B = B-max(B);
B(B<-40) = NaN;
B = B+40;


%% Linear array, direct summation without approximation, following Dobbins 2007
d = -0.002;
n = 5;
freq = 35*1e3;
c = 344;
% d = -0.01;
% n = 22;
% freq = 120*1e3;
% c = 1500;
x_pos = 0:d:(n-1)*d;
% x_pos = x_pos-mean(x_pos);
y_pos = zeros(size(x_pos));
r = 30;
theta = 0:0.05:180;
[x,y] = pol2cart(theta/180*pi,r);
k = 2*pi*freq/c;
rn = sqrt((repmat(x,length(x_pos),1)-repmat(x_pos',1,length(x))).^2 +...
     (repmat(y,length(y_pos),1)-repmat(y_pos',1,length(y))).^2);
r = sqrt((x-x_pos(1)).^2+(y-y_pos(1)).^2);
r_diff = rn-repmat(r,size(rn,1),1);
r_diff_approx = (0:n-1)'*d*cos(theta/180*pi);
b = sum(1./rn.*exp(-1j*k*r_diff),1);
B = 20*log10(abs(b));
B = B-max(B);
B(B<-40) = NaN;
B = B+40;
% mu = 0;
% sigma = abs(d);
% taper = normpdf(x_pos,mu,sigma);
% b_taper = sum(1./rn.*repmat(taper',1,size(rn,2)).*exp(-1j*k*r_diff),1);
% B_taper = 20*log10(abs(b_taper));
% B_taper = B_taper-max(B_taper);
% B_taper(B_taper<-40) = NaN;
% B_taper = B_taper+40;


%% Planar array, direct summation without approximation
d = -0.002;
n = 10;
freq = 35*1e3;
c = 344;
x_pos = 0:d:(n-1)*d;
x_pos = x_pos-mean(x_pos);
y_pos = zeros(size(x_pos));
x_pos = repmat(x_pos,1,3);
y_pos = [y_pos,y_pos+0.005,y_pos-0.005];
z_pos = zeros(size(x_pos));
r = 3;
% theta = 0:1:90;
% phi = 0:1:180;
[theta,phi] = meshgrid(0:3:90,30:3:210);
[x,y,z] = sph2cart(phi(:)/180*pi,theta(:)/180*pi,r);
k = 2*pi*freq/c;
rn = sqrt((repmat(x',length(x_pos),1)-repmat(x_pos',1,length(x))).^2 +...
     (repmat(y',length(y_pos),1)-repmat(y_pos',1,length(y))).^2 +...
     (repmat(z',length(z_pos),1)-repmat(z_pos',1,length(z))).^2);
r = sqrt((x'-x_pos(1)).^2+(y'-y_pos(1)).^2+(z'-z_pos(1)).^2);
r_diff = rn-repmat(r,size(rn,1),1);
b = sum(1./rn.*exp(-1j*k*r_diff),1);
B = 20*log10(abs(b));
B = B-max(B);
B(B<-40) = NaN;
B = B+40;

% Project lat-lon to map projection distance
[lat1,lon1] = rotatem(theta,phi,[90 90],'forward','degree');
[lat1,lon1] = rotatem(lat1,lon1,[0 -60],'forward','degree');
mstruct = defaultm('eckert4');
mstruct = defaultm(mstruct);
[xq,yq] = mfwdtran(mstruct,lat1,lon1);


%% Linear array, brute force summation
d = -0.01;
n = 21;
x_pos = 0:d:n*d;
x_pos = x_pos-mean(x_pos);
y_pos = zeros(size(x_pos));
freq = 120*1e3;
c = 1500;
theta = 0:0.5:180;
k = 2*pi*freq/c;
D = (0:n-1)*d;
b = sum(exp(-1j*k*D'*cos(theta/180*pi)),1);
B = 20*log10(abs(b));
B = B-max(B);
B(B<-40) = NaN;
B = B+40;


%% Dolphin teeth receiving array from Dobbins 2007
d = -0.01;
n = 21;
x_pos = 0:d:n*d;
% x_pos = x_pos-mean(x_pos);
y_pos = zeros(size(x_pos));
freq = 120*1e3;
c = 1500;
k = 2*pi*freq/c;
r = 30;  % in the far-field
theta = 0:0.5:180;
[x,y] = pol2cart(theta/180*pi,r);
rn = sqrt((repmat(x,length(x_pos),1)-repmat(x_pos',1,length(x))).^2+...
          (repmat(y,length(y_pos),1)-repmat(y_pos',1,length(y))).^2);
r = sqrt((x-x_pos(1)).^2+(y-y_pos(1)).^2);
fac1 = rn-repmat(r,size(rn,1),1);
b = sum(1./rn.*exp(-1j*k*fac1),1);
B = 20*log10(abs(b));
B = B-max(B);
B(B<-40) = NaN;
B = B+40;

clear
d = -0.01;
n = 22;
x_pos = 0:d:(n-1)*d;
y_pos = zeros(size(x_pos));
freq = 120*1e3;
c = 1500;
k = 2*pi*freq/c;
r = 3;  % in the far-field
theta = 0:0.5:180;
[x,y] = pol2cart(theta/180*pi,r);
rn = sqrt((repmat(x,length(x_pos),1)-repmat(x_pos',1,length(x))).^2+...
          (repmat(y,length(y_pos),1)-repmat(y_pos',1,length(y))).^2);
r = sqrt((x-x_pos(1)).^2+(y-y_pos(1)).^2);
tau = ((1:n)-1)*d/c;
omega = 2*pi*freq;
fac1 = rn-repmat(r,size(rn,1),1);
fac2 = repmat(omega*tau,size(fac1,2),1)';
b = sum(1./rn.*exp(-1j*(k*fac1-fac2)),1);
B = 20*log10(abs(b));
B = B-max(B);
B(B<-40) = NaN;
B = B+40;


%% Linear array with multiple elements
addpath('C:\Users\Wu-Jung Lee\Dropbox\0_CODE\MATLAB\array_factor')

x_pos = linspace(-0.01,0.01,11);
% y_pos = mean(diff(x_pos))*6*ones(size(x_pos));
% y_pos(1:2:end) = -y_pos(1:2:end);
y_pos = ones(size(x_pos));
w = ones(1,length(x_pos)); % 1xP vector of weighting factors
% w = normpdf(x_pos,0,0.005);

f = 50e3;
c = 344;

tt = -180:0.5:180;
% pp = 0:0.5:180;
pp = 0;
[W, theta, phi, k_x, k_y] = array_factor(x_pos, y_pos, w, f, c, tt, pp);
W = 20*log10(W);
[tt,pp] = meshgrid(theta,phi);
tt = tt/pi*180;
pp = pp/pi*180;
tt = tt-90;  % shift theta due to coord system definition                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
[tt_rot,pp_rot] = rotatem(tt,pp,[90 0],'forward','degrees');
mstruct = defaultm('ortho');
mstruct = defaultm(mstruct);
[xq_rot,yq_rot] = mfwdtran(mstruct,tt_rot,pp_rot);

figure
% contourf(xq_rot,yq_rot,W,[-30:3:0],'edgecolor','none');
surf(xq_rot,yq_rot,W','edgecolor','none');
% caxis([-30 0])
colormap(parula(9))
view([0 0 1]);
axis equal

