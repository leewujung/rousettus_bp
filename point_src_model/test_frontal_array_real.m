%% Planar array, direct summation without approximation

clear
freq_all = (20:5:50)*1e3;
% freq_all = 35e3;

c = 344;
xyz = [-2   -5    -2
       -2   5     -2 
       0    -4    1
       0    4     1
       2.5  -3    -3
       2.5  3     -3
       4    -1.7  1
       4    1.7   1
       5.5  0     -1]*1e-3;
xyz = xyz(:,[2 3 1]);
xyz = bsxfun(@minus,xyz,mean(xyz,1));
x_pos = xyz(:,1)';
y_pos = xyz(:,2)';
z_pos = xyz(:,3)';

r = 3;
[theta,phi] = meshgrid(-20:3:90,-180:3:180);
[x,y,z] = sph2cart(phi(:)/180*pi,theta(:)/180*pi,r);

w_pos = sort(repmat(normpdf(0:length(x_pos)/2,3),1,2),'ascend');
w_pos(end) = [];
w_pos = ones(1,length(x_pos));
src = [0.001,0.0015,-0.005];
steer_d = sqrt((src(1)-x_pos).^2+(src(2)-y_pos).^2+(src(3)-z_pos).^2);
steer_d = steer_d-min(steer_d);
% steer_d = zeros(size(x_pos));

% Plot
figure
plot3(x,y,z,'.');
hold on
plot3(x_pos*100,y_pos*100,z_pos*100,'-o');
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
    B = reshape(B,size(theta));
%     B(phi<30&phi>-30&theta<30&theta>-30) = 50;

    % Project lat-lon to map projection distance
    [lat1,lon1] = rotatem(theta,phi,[90 -90],'forward','degree');  % rotate the first time to go to teeth coord
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
axis equal
colormap(jet(length(freq_all)))
colorbar('Ticks',linspace(0,1,length(freq_all)),'TickLabels',{num2str(freq_all'/1e3)})


