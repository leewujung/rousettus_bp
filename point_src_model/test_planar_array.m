%% Planar array, direct summation without approximation

clear
freq_all = (20:5:50)*1e3;
% freq_all = 50e3;

d = -0.004;
n = 5;
c = 344;
x_pos = 0:d:(n-1)*d;
x_pos = x_pos-mean(x_pos);
y_pos = zeros(size(x_pos));
x_pos = [x_pos,x_pos+d/2];
y_pos = [y_pos,y_pos+0.0045];
[x_pos,isort] = sort(x_pos);
y_pos = y_pos(isort);
z_pos = zeros(size(x_pos));
w_pos = sort(repmat(normpdf(0:length(x_pos)/2-1,3),1,2));
w_pos(1) = sort(repmat(normpdf(0:length(x_pos)/2-1,3),1,2));
% w_pos = ones(1,length(x_pos));
% steer = 0;
src = [0.004,0.002,0.005];
steer_d = sqrt((src(1)-x_pos).^2+(src(2)-y_pos).^2+(src(3)-z_pos).^2);
steer_d = steer_d-min(steer_d);
r = 3;
[theta,phi] = meshgrid(0:2:180,-90:2:270);
[x,y,z] = sph2cart(phi(:)/180*pi,theta(:)/180*pi,r);

for iF = 1:length(freq_all)
    freq = freq_all(iF);
    omega = 2*pi*freq;
    k = 2*pi*freq/c;
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
%     B = B+30;

    % Project lat-lon to map projection distance
    [lat1,lon1] = rotatem(theta,phi,[90 90],'forward','degree');
    mstruct = defaultm('eckert4');
    mstruct = defaultm(mstruct);
    [xq,yq] = mfwdtran(mstruct,lat1,lon1);
    
    % Get 3dB contour
    figure
    % surf(xq,yq,reshape(B,size(xq)),'edgecolor','none')
    [C,h] = contour(xq,yq,reshape(B,size(xq)),0:-3:-30,'fill','on');
    view([0 0 1])
    colormap(parula(10))
    title(sprintf('Freq = %2.0f kHz',freq/1e3));
    
    Cout = parse_contour_output(C);
    c3db_xy = [];
    for iT=1:length(Cout)  % in case contour break into pieces
        if Cout(iT).Level == -3
            c3db_xy = [c3db_xy; Cout(iT).X',Cout(iT).Y'];
        end
    end
%     pause
%     close
    c3db_all{iF} = c3db_xy;
end


% colorset = jet(length(freq_all));
%  
% fig_3db = figure;
% hold on
% for iF=1:length(freq_all)
%     plot3(c3db_all{iF}(:,1),c3db_all{iF}(:,2),repmat(freq_all(iF)/1e3*0.01,size(c3db_all{iF},1)),...
%           'linewidth',2,'color',colorset(iF,:));
% end
% axis equal
% colormap(jet(length(freq_all)))
% colorbar('Ticks',linspace(0,1,length(freq_all)),'TickLabels',{num2str(freq_all'/1e3)})


