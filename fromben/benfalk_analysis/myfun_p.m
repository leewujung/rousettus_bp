function F = myfun_p(x,speaker_pos,delay_mic,c)
% delay_mic  [1 x #speaker]

fac = sum((repmat(x,size(speaker_pos,1),1)-speaker_pos).^2,2);

% F = sqrt((c*delay_mic').^2-fac);

if size(delay_mic,2)~=1
    delay_mic = delay_mic';
end

F = (c*delay_mic).^2-fac;