
theta=-pi/2:pi/200:pi/2;
phi=0:pi/500:2*pi;

[TT,PP] = meshgrid(theta,phi);
[XX,YY,ZZ] = sph2cart(reshape(TT,1,[]),reshape(PP,1,[]),ones(1,numel(TT)));

yvec = [0,1,0];
pol_angle = acos([XX;YY;ZZ]'*yvec');
az = atan(-ZZ./XX);


c = 344;
freq = 55e3;
k = 2*pi*freq/c;
a = 4e-3;
rr = 2*besselj(1,k*a*sin(pol_angle))./(k*a*sin(pol_angle));
rr(pol_angle==0) = 1;

rr2 = 2*besselj(1,k*a*sin(az))./(k*a*sin(az));
rr2(az==0) = 1;
% rr = rr.*rr2.';

[XXbeam,YYbeam,ZZbeam] = sph2cart(reshape(TT,1,[])',reshape(PP,1,[])',abs(rr));
% [XXbeam,YYbeam,ZZbeam] = sph2cart(reshape(TT,1,[])',reshape(PP,1,[])',abs(rr2)');
% [XXbeam,YYbeam,ZZbeam] = sph2cart(reshape(TT,1,[])',reshape(PP,1,[])',20*log10(abs(rr))-min(20*log10(abs(rr))));

% figure;
% plot3(XXbeam,YYbeam,ZZbeam,'b.')
% grid on
% axis equal

XXbeam = reshape(XXbeam,size(TT,1),[]);
YYbeam = reshape(YYbeam,size(TT,1),[]);
ZZbeam = reshape(ZZbeam,size(TT,1),[]);

figure;
surf(XXbeam,YYbeam,ZZbeam,ones(size(ZZbeam)));
grid on
% axis equal
