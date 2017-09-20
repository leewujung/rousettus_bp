% 2015 11 29  Test Fourier transform technique for bp calculation

d = -0.002;
n = 6;
freq = 35*1e3;
c = 344;
x_pos = 0:d:(n-1)*d;
x_pos = x_pos-mean(x_pos);
y_pos = zeros(size(x_pos));
r = 3;
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

c = 344;  % sound speed [m/s]
freq = 35e3;  % frequency [Hz]
lambda = c/freq;  % wavelength [m]
T = 1e-3;  % unit length: 1 mm
n = 6;  % 5 elementbs
d = 0.002;  % [m]
L = n*d;
w = ones(round(L/T),1);
w = [w;zeros(1000,1)];  % zero-padding
W = ifftshift(fft(w));
W = 20*log10(abs(W));
W = W-max(W);
W(W<-40) = NaN;
W = W+40;
omega = linspace(0,1/T/2,length(W)/2);
omega = [-fliplr(omega),omega];
beta = lambda*omega;
theta = asin(beta);
m = -10:10;
delta_beta = m*lambda/d;

H = conv(W,)
