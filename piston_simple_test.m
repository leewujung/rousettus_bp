
c = 343;
freq = (35:1:90)*1e3;
bp_info.k = 2*pi*freq'/c;
pol_angle = 0:pi/180:pi/2;
bp_info.a = 4e-3;

mic_dB = 20*log10(abs(2*besselj(1,bp_info.k*bp_info.a*sin(pol_angle))./...
                     (bp_info.k*bp_info.a*sin(pol_angle))));
                 
imagesc(pol_angle/pi*180,freq/1e3,mic_dB);
xlabel('Angle (deg)');
ylabel('Frequency (kHz)');
title(sprintf('Aperture radius = %3.1f mm',bp_info.a*1e3));
colorbar
caxis([-30 0])