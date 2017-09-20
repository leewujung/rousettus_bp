%% Planar array, direct summation without approximation

clear
freq_all = (20:5:50)*1e3;
% freq_all = 35e3;

angle = 14.5;  % jaw angle to midline [deg]
x_pos1 = [3, 5.5, 9, 11, 15,...
         2, 4.5, 7.5, 10, 14]*1e-3;
y_pos1 = tan(angle/180*pi)*x_pos1;
y_pos1 = y_pos1-y_pos1(end);
z_pos1 = [0, 0, 0, 0, 0,...
         -3.5, -4, -5, -6, -6]*1e-3;
% x_pos2 = x_pos1([3,4,8,9]);
% y_pos2 = -y_pos1([3,4,8,9]);
% z_pos2 = z_pos1([3,4,8,9]);
x_pos2 = x_pos1([1:4,6:9]);
y_pos2 = -y_pos1([1:4,6:9]);
z_pos2 = z_pos1([1:4,6:9]);
xyz(:,1) = [x_pos1,x_pos2];
xyz(:,2) = [y_pos1,y_pos2];
xyz(:,3) = [z_pos1,z_pos2];
% xyz(:,1) = x_pos1;
% xyz(:,2) = y_pos1;
% xyz(:,3) = z_pos1;

c = 344;
xyz = xyz(:,[2 3 1]);
xyz = sortrows(xyz,3);
xyz = bsxfun(@minus,xyz,mean(xyz,1));
xxx_pos = xyz(:,1)';
yyy_pos = xyz(:,2)';
zzz_pos = xyz(:,3)';

r = 3;
[theta,phi] = meshgrid(-20:3:90,-180:3:180);
[x,y,z] = sph2cart(phi(:)/180*pi,theta(:)/180*pi,r);

% w_pos = sort(repmat(normpdf(0:length(x_pos)/2,3),1,2),'ascend');
% w_pos(end) = [];
w_pos = ones(1,length(xxx_pos));
src = [0.002,0,-0.01];
steer_ddd = sqrt((src(1)-xxx_pos).^2+(src(2)-yyy_pos).^2+(src(3)-zzz_pos).^2);
steer_ddd = steer_ddd-min(steer_ddd);
% steer_ddd = zeros(size(xxx_pos));
% steer_ddd = sort(repmat((0:length(xxx_pos)/2)*0.5,1,2),'ascend');
% steer_ddd(end) = [];

% Plot
figure
plot3(x,y,z,'.','markersize',5);
hold on
plot3(xxx_pos*100,yyy_pos*100,zzz_pos*100,'-o');
plot3(src(1)*100,src(2)*100,src(3)*100,'r*')
quiver3(0,0,0,2,0,0,'k','linewidth',2)
quiver3(0,0,0,0,2,0,'k','linewidth',2)
quiver3(0,0,0,0,0,2,'k','linewidth',2)
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on
axis equal


fig_geo = figure;
fig_con = figure;
for iF = 1:length(freq_all)
    freq = freq_all(iF);
    switch iF
        case {1,2}
            idx1 = 9:length(xxx_pos);
        case {3,4}
            idx1 = 5:length(xxx_pos);
        case {5,6,7}
            idx1 = 1:length(xxx_pos);
    end
    idx2 = ones(1,length(xxx_pos));
    if src(1)>0
        idx2([1,3,5,7,9,11])=0;
    else
        idx2([2,4,6,8,10,12])=0;
    end
    idx = intersect(idx1,find(idx2));
%     idx = 1:length(xxx_pos);
    x_pos = xxx_pos(idx);
    y_pos = yyy_pos(idx);
    z_pos = zzz_pos(idx);
    steer_d = steer_ddd(idx);
    w_pos = ones(1,length(x_pos));
    omega = 2*pi*freq;
    k = 2*pi*freq/c;
    rn = sqrt((repmat(x',length(x_pos),1)-repmat(x_pos',1,length(x))).^2 +...
        (repmat(y',length(y_pos),1)-repmat(y_pos',1,length(y))).^2 +...
        (repmat(z',length(z_pos),1)-repmat(z_pos',1,length(z))).^2);
    r = sqrt((x'-x_pos(1)).^2+(y'-y_pos(1)).^2+(z'-z_pos(1)).^2);
    r_diff = rn-repmat(r,size(rn,1),1);
    fac1 = 1./rn;
%     fac2 = exp( -1j*(k*r_diff) );
    fac2 = exp( -1j*(k*(r_diff+repmat(steer_d(:),1,size(r_diff,2)))) );
    fac3 = repmat(w_pos(:),1,size(r_diff,2));
    b = sum(fac1.*fac2.*fac3,1);
    B = 20*log10(abs(b));
    B = B-max(B);
    B(B<-30) = NaN;
    B = reshape(B,size(theta));
%     B(phi<30&phi>-30&theta<30&theta>-30) = 50;

    % Project lat-lon to map projection distance
    [lat1,lon1] = rotatem(theta,phi,[90 -90],'forward','degree');  % rotate the first time to go to teeth coord
    lon1 = -lon1;
%     lon1 = 90-lon1;  % flip +phi direction and rotate to bat coord
%     lon1 = lon1-14.5;  % correct for teeth angle to bat coord
    mstruct = defaultm('eckert4');
    mstruct = defaultm(mstruct);
    [xq,yq] = mfwdtran(mstruct,lat1,lon1);
    
    % Geoshow
    figure(fig_geo);
    subplot(2,4,iF);
    axesm eckert4;
    framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-120 120]);
    gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
	geoshow(lat1,lon1,B,'displaytype','texturemap');
    caxis([-30 0])

    % Get 3dB contour
    figure(fig_con);
    subplot(2,4,iF);
    [C,h] = contour(xq,yq,B,0:-3:-30,'fill','on');
    view([0 0 1])
    colormap(parula(10))
    title(sprintf('Freq = %2.0f kHz',freq/1e3));
    axis equal
    
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


colorset = jet(length(freq_all));
 
fig_3db = figure;
hold on
for iF=1:length(freq_all)
    plot3(c3db_all{iF}(:,1),c3db_all{iF}(:,2),repmat(freq_all(iF)/1e3*0.01,size(c3db_all{iF},1)),...
          'linewidth',2,'color',colorset(iF,:));
end
colormap(jet(length(freq_all)))
colorbar('Ticks',linspace(0,1,length(freq_all)),'TickLabels',{num2str(freq_all'/1e3)})
axis equal


