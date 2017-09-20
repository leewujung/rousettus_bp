function mic_dB = model_beam_dB(bp_info,beam_aim,mic_azel)
% Predicted beam pattern based on piston or gaussian model
% INPUT
%   bp_info    prescribed bp info
%   beam_aim   beam aim from experimental data
%   mic_azel   az,el position of mic from experimental data [rad]
% OUTPUT
%   mic_dB     predicted mic receiving intensity
%
% Wu-Jung Lee | leewujung@gmail.com

[aim_xyz(:,1),aim_xyz(:,2),aim_xyz(:,3)] = sph2cart(beam_aim(:,1),beam_aim(:,2),1);
[mic_xyz(:,1),mic_xyz(:,2),mic_xyz(:,3)] = sph2cart(mic_azel(:,1),mic_azel(:,2),1);

pol_angle = acos(aim_xyz*mic_xyz');
switch bp_info.type
    case 'piston'
%         pol_angle(pol_angle>pi/2) = NaN;
        mic_dB = 20*log10(abs(2*besselj(1,bp_info.k*bp_info.a*sin(pol_angle))./...
                     (bp_info.k*bp_info.a*sin(pol_angle))));
        mic_dB(pol_angle==0) = 0;
    case 'gaussian'
        mic_dB = 20*log10(normpdf(pol_angle,bp_info.mu,bp_info.sigma));
end

mic_dB(abs(pol_angle)>pi/2) = ...   % "flatten" out the 2nd mainlobe arising from wrapping around piston model
    20*log10(abs(2*besselj(1,bp_info.k*bp_info.a*sin(pi/2))./...
    (bp_info.k*bp_info.a*sin(pi/2))));
% mic_dB(isnan(mic_dB)) = -Inf;



