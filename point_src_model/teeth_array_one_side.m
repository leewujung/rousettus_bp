%% Planar array, direct summation without approximation
clear
d = -0.002;
n = 10;
% freq = 35*1e3;
freq_all = (20:5:50)*1e3;
c = 344;

x_pos = [3, 5.5, 9, 11, 15,...
         2, 4.5, 7.5, 10, 14]*1e-3*1.2;
z_pos = [0, 0, 0, 0, 0,...
         -3.5, -4, -5, -6, -6]*1e-3;
x_pos = x_pos-mean(x_pos);
z_pos = z_pos-mean(z_pos);
y_pos = zeros(size(x_pos));
     
% x_pos = 0:d:(n-1)*d;
% x_pos = x_pos-mean(x_pos);
% z_pos = zeros(size(x_pos));
% x_pos = repmat(x_pos,1,3);
% z_pos = [z_pos,z_pos+0.002,z_pos-0.002];
% y_pos = zeros(size(x_pos));

r = 3;
[theta,phi] = meshgrid(-90:3:90,-100:3:180);
[x,y,z] = sph2cart(phi(:)/180*pi,theta(:)/180*pi,r);

x_mark = x(phi<30&phi>-30&theta<30&theta>-30);
y_mark = y(phi<30&phi>-30&theta<30&theta>-30);
z_mark = z(phi<30&phi>-30&theta<30&theta>-30);

src = [-.01,-0.005,0];
steer_d = sqrt((src(1)-x_pos).^2+(src(2)-y_pos).^2+(src(3)-z_pos).^2);
steer_d = steer_d-min(steer_d);

w_pos = ones(1,length(x_pos));


% Plot
figure
plot3(x,y,z,'.');
hold on
plot3(x_pos*100,y_pos*100,z_pos*100,'-o');
plot3(x_mark,y_mark,z_mark,'o');
plot3(src(1)*100,src(2)*100,src(3)*100,'r*');
quiver3(0,0,0,2,0,0,'k','linewidth',2)
quiver3(0,0,0,0,2,0,'k','linewidth',2)
quiver3(0,0,0,0,0,2,'k','linewidth',2)
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on
axis equal


% Calculate beampattern

fig_geo = figure;
fig_con = figure;
for iF = 1:length(freq_all)
    freq = freq_all(iF);

k = 2*pi*freq/c;
% rn = sqrt((repmat(x',length(x_pos),1)-repmat(x_pos',1,length(x))).^2 +...
%           (repmat(y',length(y_pos),1)-repmat(y_pos',1,length(y))).^2 +...
%           (repmat(z',length(z_pos),1)-repmat(z_pos',1,length(z))).^2);
% r = sqrt((x'-x_pos(1)).^2+(y'-y_pos(1)).^2+(z'-z_pos(1)).^2);
% r_diff = rn-repmat(r,size(rn,1),1);
% b = sum(1./rn.*exp(-1j*k*r_diff),1);

rn = sqrt((repmat(x',length(x_pos),1)-repmat(x_pos',1,length(x))).^2 +...
          (repmat(y',length(y_pos),1)-repmat(y_pos',1,length(y))).^2 +...
          (repmat(z',length(z_pos),1)-repmat(z_pos',1,length(z))).^2);
r = sqrt((x'-x_pos(1)).^2+(y'-y_pos(1)).^2+(z'-z_pos(1)).^2);
r_diff = rn-repmat(r,size(rn,1),1);
fac1 = 1./rn;
fac2 = exp( -1j*(k*(r_diff+repmat(steer_d(:),1,size(r_diff,2)))) );
fac3 = repmat(w_pos(:),1,size(r_diff,2));
b = sum(fac1.*fac2.*fac3,1);

B = 20*log10(abs(b));
B = B-max(B);
B(B<-30) = NaN;

B = reshape(B,size(theta));
% B(phi<30&phi>-30&theta<30&theta>-30) = 30;


% Projection
% Note: "-phi" is used instead of phi to flip the left/right orientation
figure;
axesm eckert4;
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
geoshow(theta,-phi,B,'displaytype','texturemap');
caxis([-30 0])

