%Piston model
%           |2*J1(k*a*sin(theta))|
% RP(theta)=|--------------------|
%           |   k*a*sin(theta)   |
%RP(theta) = ratio between on-axis pressure and pressure at a given angle theta
%J1 = first-order Bessel function
%k = 2pi/lambda
%lambda = wavelength (c/f)
%c = speed of sound in air (~343 m/s)
%a = radius of the sound emitter
%http://www.acs.psu.edu/drussell/Demos/BaffledPiston/BaffledPiston.html


close all;
clear;

save_video=0;

plot_contour_3D=1;     map_proj=1;

plot_3D=0;
plot_3D_kaushik=0;
plot_surface=0; %alternate is a scatter
plot_3D_slices=0;
plot_2D=1;
plot_atm_atten=0;

angle_step_size=pi/1e3;
c=344; %speed of sound
th=-pi/2:angle_step_size:pi/2;
phi=0:angle_step_size:2*pi;
freq=25e3:.5e3:90e3;
ah=5.9/1e3; %mm to meters - aperture, use the same for both vertical and horiz
av=5.9/1e3;

dblim=30;
colordef black
% colordef white

if plot_contour_3D
  figure(1); clf;
  set(gcf, 'pos',[180 60 600 600],'color','k')
    
  if save_video
    writerObj = VideoWriter('F:\piston-model-contour-3D','Uncompressed AVI');
    writerObj.FrameRate=20;
    %   writerObj.Quality=90;
    open(writerObj);
  end
  
  if map_proj
    axesm ortho;
    framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-120 120]);
    gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
    axis off
    tightmap
    cb=colorbar('location','southoutside');
    set(cb,'visible','off')
  end
  
  for f=1:length(freq)
    
%           |2*J1(k*a*sin(theta))|
% RP(theta)=|--------------------|
%           |   k*a*sin(theta)   |
    
    lambda=c/freq(f);
    k=2*pi/lambda;
    kasinh=k.*ah.*sin(th);
    Jh = besselj(1,kasinh);
    Rh= abs(2.*Jh) ./ abs(kasinh);
    Rh(th==0)=1;
    
    [x,y,z] = pol2cart(th,Rh,zeros(size(th)));
    xyzrot={};
    for pp=1:5:length(phi)
      rotx=[1 0 0;0 cos(phi(pp)) -sin(phi(pp));0 sin(phi(pp)) cos(phi(pp))];
      xyzrot{end+1}=(rotx*[x;y;z]);
      %       scatter3(xyzrot{pp}(:,1),xyzrot{pp}(:,2),xyzrot{pp}(:,3),'filled')
      %       hold on;
    end
    
    xyz=cell2mat(xyzrot)';
    x=xyz(:,1);
    y=xyz(:,2);
    z=xyz(:,3);
    
    xresh=reshape(x,length(th),[]);
    yresh=reshape(y,length(th),[]);
    zresh=reshape(z,length(th),[]);
    [az,el,vq]=cart2sph(xresh,yresh,zresh);
    
    vqnorm=20*log10(vq);
    
    cla
    if map_proj
      axesm ortho;
      framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-120 120]);
      gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
%       axis off
      geoshow(el/pi*180,az/pi*180,vqnorm,'displaytype','texturemap');
      contourm(el/pi*180,az/pi*180,vqnorm,-3,'w','linewidth',2);
      
    else
      contourf(az,el,vqnorm,'LevelStep',3,'LineStyle','none')
      axis equal;
    end
    caxis([-30 0])
    colorbar('southoutside');
    
    title([num2str(round(freq(f)/1e3)) ' kHz ' ...
      char(10) ' @ ' num2str(ah*1e3) ' mm aperture'])
    
    if save_video
      writeVideo(writerObj,getframe(gcf));
    else
      drawnow;
    end
    
  end
  
  if save_video
    close(writerObj);
    fname=[writerObj.Path writerObj.Filename];
    compress_video(fname,1,1,'.\');
  end
  
end


if plot_3D
  figure(1); clf;
  set(gcf, 'pos',[180 60 600 740],'color','k')
  
  %3D polar
  if save_video
    writerObj = VideoWriter('F:\piston-model-3D','Uncompressed AVI');
    writerObj.FrameRate=20;
    %   writerObj.Quality=90;
    open(writerObj);
  end
  for f=1:length(freq)
    lambda=c/freq(f);
    k=2*pi/lambda;
    
    kasinh=k.*ah.*sin(th);
    Jh = besselj(1,kasinh);
    Rh= abs(2.*Jh) ./ abs(kasinh);
    Rh(th==0)=1;
    
    [x,y,z] = pol2cart(th,Rh,zeros(size(th)));
    %     scatter3(x,y,z)
    %phi is a rotation around the axis of the beam
    xyzrot={};
    for pp=1:5:length(phi)
      rotx=[1 0 0;0 cos(phi(pp)) -sin(phi(pp));0 sin(phi(pp)) cos(phi(pp))];
      xyzrot{end+1}=(rotx*[x;y;z]);
      %       scatter3(xyzrot{pp}(:,1),xyzrot{pp}(:,2),xyzrot{pp}(:,3),'filled')
      %       hold on;
    end
    
    xyz=cell2mat(xyzrot)';
    x=xyz(:,1);
    y=xyz(:,2);
    z=xyz(:,3);
    
    if plot_surface
      xresh=reshape(x,length(th),[]);
      yresh=reshape(y,length(th),[]);
      zresh=reshape(z,length(th),[]);
      surf(xresh,yresh,zresh);
      if f==1
        colormap gray
      end
      light('Position',[50 -15 29])
      lighting phong
      shading interp
      view(-20,20)
      axis equal
      axis([0 1 -1 1 -1 1])
    else
      if freq(f)==35e3
        col='r';
      else
        col='w';
      end
      scatter3(x,y,z,3,col,'filled')
      view(-20,20)
      axis equal
      axis([0 1 -1 1 -1 1])
    end
    
    title([num2str(round(freq(f)/1e3)) ' kHz ' ...
      char(10) ' @ ' num2str(ah*1e3) ' mm aperture'])
    %   axis([0 dblim -dblim dblim -dblim dblim])
    
    
    if save_video
      writeVideo(writerObj,getframe(gcf));
    else
      drawnow;
    end
    
    if (freq(f)==35e3 || f==length(freq))
      [vv_az,vv_el]=view;
      for vv=[vv_az:-1:vv_az-60 vv_az-60:vv_az+60 vv_az+60:-1:vv_az]
        view(vv,vv_el);
        if save_video
          writeVideo(writerObj,getframe(gcf));
        else
          drawnow;
        end
      end
    end
  end
  
  if save_video
    close(writerObj);
    compress_video('F:\piston-model-3D.avi',1,0)
  end
end



if plot_3D_slices
  figure(4); clf;
  set(gcf, 'pos',[180 60 600 740],'color','k')
  
  if save_video
    writerObj = VideoWriter('F:\piston-model-3D','Uncompressed AVI');
    writerObj.FrameRate=20;
    %   writerObj.Quality=90;
    open(writerObj);
  end
  for f=1:length(freq)
    lambda=c/freq(f);
    k=2*pi/lambda;
    kasinh=k.*ah.*sin(th);
    kasinv=k.*av.*sin(phi);
    
    Jh = besselj(1,kasinh);
    Jv = besselj(1,kasinv);
    
    Rh=abs(2.*Jh) ./ abs(kasinh);
    Rv=abs(2.*Jv) ./ abs(kasinv);
    
    R=Rh'*Rv;
    %     polar3d(R,min(th),max(th),0,1,1,'surf','cubic');
    
    for pp=1:length(phi)
      [x,y,z] = sph2cart(th,repmat(phi(pp),size(th)),R(:,pp)');
      subplot(2,1,1)
      scatter3(x,y,z,3,'w','filled')
      title(['phi: ' num2str(phi(pp)*180/pi) ' '...
        char(10) num2str(round(freq(f)/1e3)) ' kHz ' ...
        ' @ ' num2str(ah*1e3) ' mm aperture'])
      %       view(-20,20)
      axis equal
      axis([0 1 -.8 .8 ])
      view(2)
      drawnow
      
      subplot(2,1,2)
      scatter3(x,y,z,3,'w','filled')
      title(['phi: ' num2str(phi(pp)*180/pi) ' '...
        char(10) num2str(round(freq(f)/1e3)) ' kHz ' ...
        ' @ ' num2str(ah*1e3) ' mm aperture'])
      %       view(-20,20)
      axis equal
      axis([0 1 -.8 .8 -.8 .8])
      view(0,0)
      drawnow
    end
    
    if save_video
      writeVideo(writerObj,getframe(gcf));
    else
      drawnow;
    end
    
    if freq(f)==35e3 || f==length(freq)
      [vv_az,vv_el]=view;
      for vv=[vv_az:-1:vv_az-60 vv_az-60:vv_az+90 vv_az+90:-1:vv_az]
        view(vv,vv_el);
        if save_video
          writeVideo(writerObj,getframe(gcf));
        else
          drawnow;
        end
      end
    end
  end
  
  if save_video
    close(writerObj);
    compress_video('F:\piston-model-3D.avi',1,0)
  end
end



if plot_2D
  figure(2); clf;
  set(gcf, 'pos',[60 60 1080 550],'color','k')
  
  if save_video
    writerObj = VideoWriter('piston-model','MPEG-4');
    writerObj.FrameRate = 10;
    writerObj.Quality=90;
    open(writerObj);
  end
  for f=1:length(freq)
    lambda=c/freq(f);
    k=2*pi/lambda;
    J = besselj(1,k*ah*sin(th));
    R=abs(2*J) ./ abs(k*ah*sin(th));
    
    Rdb=mag2db(R);
    Rdb(Rdb<-dblim)=nan;
    h=polar(th,dblim+Rdb);
    if freq(f)==35e3
      col='r';
    else
      col='y';
    end
    set(h,'linewidth',3,'color',col);
    xlim([0 dblim])
    ylim([-dblim dblim])
    view(-90,90);
    
    aaa=findall(gca,'type','text');
    rtick={'  30','  20','  10'};
    for rt=1:length(rtick)
      indx=~cellfun('isempty',strfind({aaa.String},rtick{rt}));
      set(aaa(indx),'String',[num2str(str2double(rtick{rt})-dblim) ' dB'],'fontsize',14)
    end
    
    aaa=findall(gca,'type','text');
    rtick={'0','30','60','90'};
    for rt=1:length(rtick)
      indx=strcmp({aaa.String},rtick{rt});
      set(aaa(indx),'String',num2str(str2double(rtick{rt})*-1),'fontsize',14)
    end
    
    aaa=findall(gca,'type','text');
    rtick={'330','300','270'};
    for rt=1:length(rtick)
      indx=strcmp({aaa.String},rtick{rt});
      set(aaa(indx),'String',num2str(360-str2double(rtick{rt})),'fontsize',14)
    end
    
    ylabel([num2str(round(freq(f)/1e3)) ' kHz ' ...
      char(10) ' @ ' num2str(ah*1e3) ' mm aperture'],'fontsize',14)
    if save_video
      if freq(f)==35e3
        for kk=1:3*writerObj.FrameRate-1
          writeVideo(writerObj,getframe(gcf));
        end
      else
        writeVideo(writerObj,getframe(gcf));
      end
    else
      drawnow;
      %   pause(.05)
    end
  end
  
  if save_video
    close(writerObj);
  end
  
end



%kaushik version
if plot_3D_kaushik
  
  theta=-pi/2:pi/100:pi/2;
  phi=0:pi/100:2*pi;
  a = 5.9e-3;
  freq=20e3:.5e3:90e3;
  
  [TT,PP] = meshgrid(theta,phi);
  [XX,YY,ZZ] = sph2cart(reshape(TT,1,[]),reshape(PP,1,[]),ones(1,numel(TT)));
  %   az = atan(-ZZ./XX);
  
  yvec = [0,1,0];
  pol_angle = acos([XX;YY;ZZ]'*yvec');
  
  yvec = [0,1,1];
  pol_angle1 = acos([XX;YY;ZZ]'*yvec'/norm(yvec));
  
  yvec = [0,1,-1];
  pol_angle2 = acos([XX;YY;ZZ]'*yvec'/norm(yvec));
  
  
  
  figure(1); clf;
  set(gcf, 'pos',[180 60 600 740],'color','k')
  
  %3D polar
  if save_video
    writerObj = VideoWriter('F:\piston-model-3D','Uncompressed AVI');
    writerObj.FrameRate=20;
    %   writerObj.Quality=90;
    open(writerObj);
  end
  for f=1:length(freq)
    k = 2*pi*freq(f)/c;
    
    %     rr = 2*besselj(1,k*a*sin(pol_angle))./(k*a*sin(pol_angle));
    %     rr(pol_angle==0) = 1;
    
    rr1 = 2*besselj(1,k*a*sin(pol_angle1))./(k*a*sin(pol_angle1));
    rr1(pol_angle1==0) = 1;
    
    rr2 = 2*besselj(1,k*a*sin(pol_angle2))./(k*a*sin(pol_angle2));
    rr2(pol_angle2==0) = 1;
    
    rr=rr1+rr2;
    
    %     rr2 = 2*besselj(1,k*a*sin(az))./(k*a*sin(az));
    %     rr2(az==0) = 1;
    %     rr = rr.*rr2.';
    
    [XXbeam,YYbeam,ZZbeam] = sph2cart(reshape(TT,1,[])',reshape(PP,1,[])',abs(rr));
    %     [XXbeam,YYbeam,ZZbeam] = sph2cart(reshape(TT,1,[])',reshape(PP,1,[])',20*log10(abs(rr))-min(20*log10(abs(rr))));
    
    XXbeam = reshape(XXbeam,size(TT,1),[]);
    YYbeam = reshape(YYbeam,size(TT,1),[]);
    ZZbeam = reshape(ZZbeam,size(TT,1),[]);
    
    if plot_surface
      surf(XXbeam,YYbeam,ZZbeam,ones(size(ZZbeam)));
      if f==1
        colormap gray
      end
      light('Position',[50 -15 29])
      lighting phong
      shading interp
    else
      if freq(f)==35e3
        col='r';
      else
        col='w';
      end
      scatter3(x,y,z,3,col,'filled')
    end
    title([num2str(round(freq(f)/1e3)) ' kHz ' ...
      char(10) ' @ ' num2str(ah*1e3) ' mm aperture'])
    view(70,15)
    axis equal
    axis([-1.7 1.7 0 1 -1.7 1.7])
    %     axis([-2 2 0 2 -2 2])
    
    if save_video
      writeVideo(writerObj,getframe(gcf));
    else
      drawnow;
    end
    
    if freq(f)==35e3 || f==length(freq)
      [vv_az,vv_el]=view;
      for vv=[vv_az:-1:vv_az-60 vv_az-60:vv_az+60 vv_az+60:-1:vv_az]
        view(vv,vv_el);
        if save_video
          writeVideo(writerObj,getframe(gcf));
        else
          drawnow;
        end
      end
    end
  end
  
  if save_video
    close(writerObj);
    compress_video('F:\piston-model-3D.avi',1,0)
  end
end



%     lambda=c/freq(f);
%     k=2*pi/lambda;
%
%     delta=pi/4;
%
%     kasinh=k.*ah.*sin(th);
%     Jh = besselj(1,kasinh);
%     Rh=(2.*Jh./kasinh);
%     Rh(th==0)=1;
%
%     kasinvent=k.*av.*sin(phi-delta);
%     Jvent = besselj(1,kasinvent);
%     Rvent=(2.*Jvent./kasinvent);
%     Rvent(phi-delta==0)=1;
%
%     kasindors=k.*av.*sin(phi+delta);
%     Jdors = besselj(1,kasindors);
%     Rdors=(2.*Jdors./kasindors);
%     Rdors(phi+delta==0)=1;
%
% %     Rv=ones(size(Rh));
%
%     Rdors=(Rh') * Rdors;
%     Rvent=(Rh') * Rvent;
%
%     R=Rdors.^2+Rvent.^2;
%     Rexpanded=R(:);
%     %   RexpdB=mag2db(Rexpanded);
%     %   RexpdB(RexpdB<-dblim)=nan;
%     thexpanded=repmat(th,size(th))';
%     phiexpanded=repmat(phi',size(phi))';
%     phiexpanded=phiexpanded(:);
%
%
%     [x,y,z] = sph2cart(thexpanded,phiexpanded,Rexpanded);
% %     phireal=(dot([x1 y1 z1],repmat([1 0 0],size(x1)),2));
% %     [x,y,z] = sph2cart(thexpanded,phireal,Rexpanded);
%
%     if plot_surface
%       xresh=reshape(x,length(th),[]);
%       yresh=reshape(y,length(th),[]);
%       zresh=reshape(z,length(th),[]);
%       surf(xresh,yresh,zresh);
%       if f==1
%         colormap gray
%       end
%       light('Position',[50 -15 29])
%       lighting phong
%       shading interp
%     else
%       if freq(f)==35e3
%         col='r';
%       else
%         col='w';
%       end
%       scatter3(x,y,z,3,col,'filled')
%     end
%
%     title([num2str(round(freq(f)/1e3)) ' kHz ' ...
%       char(10) ' @ ' num2str(ah*1e3) ' mm aperture'])
%     view(-20,20)
%     axis equal
%     axis([0 1 -1 1 -1 1])
%     %   axis([0 dblim -dblim dblim -dblim dblim])
%
%
%     if make_anim
%       writeVideo(writerObj,getframe(gcf));
%     else
%       drawnow;
%     end
%
%     if freq(f)==35e3 || f==length(freq)
%       [vv_az,vv_el]=view;
%       for vv=[vv_az:-1:vv_az-60 vv_az-60:vv_az+60 vv_az+60:-1:vv_az]
%         view(vv,vv_el);
%         if make_anim
%           writeVideo(writerObj,getframe(gcf));
%         else
%           drawnow;
%         end
%       end
%     end
%   end
%
%   if make_anim
%     close(writerObj);
%     compress_video('F:\piston-model-3D.avi',1,0)
%   end
% end


if plot_atm_atten
  %from http://www.kayelaby.npl.co.uk/general_physics/2_4/2_4_1.html
  atm_att=[10,190,280,240,190,160,130,120,100,95;12.5000000000000,210,360,340,280,240,200,180,160,140;16,230,430,470,420,360,320,280,250,230;20,260,510,600,580,520,470,420,380,350;25,300,580,740,770,730,680,620,570,520;31.5000000000000,360,670,890,990,1000,960,900,840,790;40,460,780,1100,1200,1300,1300,1300,1200,1200;50,600,940,1300,1500,1700,1700,1700,1700,1700;63,840,1200,1500,1800,2100,2200,2300,2300,2300;80,1200,1600,2000,2300,2600,2800,3000,3100,3100;100,1800,2200,2500,2900,3300,3600,3800,4000,4100];
  
  figure(3);
  semilogx(atm_att(:,1)*1e3,atm_att(:,2:end)/1e3)
  legend({'10'	'20'	'30'	'40'	'50'	'60'	'70'	'80'	'90'})
  xlabel('Frequency(Hz)')
  ylabel('Alpha (dB/m)')
  title('Atm. att. for rel. humidity')
  set(gca,'fontsize',20)
end
%formula:
%http://www.sengpielaudio.com/AirdampingFormula.htm
% There are 4 parameters for calculating air damping:
%  1) Frequency f in Hz.
%  2) Temperature T in K (Kelvin). Can be fixed to a standard room temperature of 293.15 K
%  = 20°C = 68°F.
%  3) Relative Humidity h in %. Can be fixed to a standard of 50 to 60% or whatever represents a standard RH for one's circumstances.
%  4) Atmospheric pressure p in kPa (kilopascal). Can be fixed to a standard of 101.325 kPa (standard pressure at sea level).
%
%  as = a · s [dB] total absorption at distance s
%
%  pt = pi · exp(?x · as) [Pa]
%  x = 1 / (10 · log ((exp(1))2) = ca. 0.1151
%  Delta Lt = 20 · log (pi / pt) = as [dB]
%  a = 8.686 · f2 · ((1.84 · 10?11 · (pa / pr)?1 · (T / To)1/2) + y) [dB/m]
%  y = (T / To)?5/2· (0.01275 · exp (?2239.1 / T) · (frO + f2 / frO)?1 + z)
%  z = 0.1068 · exp (?3352 / T) · (frN + f2 / frN)?1
%  frO = (pa / pr) · (24 + 4.04 · 104 · h · ((0.02 + h) / (0.391 + h)))
%  frN = (pa / pr) · (T / To)?1/2 · (9 + 280 · h · exp (?4.170 · ((T / To)?1/3?1)))
%  h = hr · ((psat / pr) / (pa / pr)) = hr · (psat / pa)
%  psat = pr · 10(?6.8346 · (To1 / T)^1.261 + 4.6151)
%
%  a ........ pure-tone sound attenuation coefficient, in dB/m, for atmospheric absorption
%  s ........ distance in m through which the sounds propagates
%  pi ....... initial sound pressure amplitude in Pa
%  pt ....... sound pressure amplitude in Pa
%  pa ...... ambient atmospheric pressure in kPa
%  pr ....... reference ambient atmospheric pressure: 101.325 kPa
%  psat ... saturation vapor pressure equals:
%  ................ International Meteorological Tables WMO-No.188 TP94
%  ................ World Meteorological Organization - Geneva Switzerland
%  T ........ ambient atmospheric temperature in K (Kelvin).
%  ........... K = 273.15 + Temperature in °C (Celsius)
%  To ...... reference temperature in K: 293.15 K (20°C)
%  To1..... triple-point isotherm temp: 273.16 K = 273.15 + 0.01 K (0.01°C)
%  h ........ molar concentration of water vapor, as a percentage
%  hr........ relative humidity as a percentage
%  f ......... frequency
%  frO ..... oxygen relaxation frequency
%  frN ..... nitrogen relaxation frequency
%  x ........ a help factor to shorten formula ? improvement by E. Desart
%  y ........ a help factor to shorten formula
%  z ........ a help factor to shorten formula



%alternative ways to plot:

%   [X,Y]=pol2cart(th+pi/2,R);
%   plot(X,Y,'linewidth',2);
%   axis equal
%   grid on
%   axis([-1 1 0 1])

%   mmpolar(th,mag2db(R),...
%     'style','compass','Rlimit',[-30 0],'tlimit',[th(1) th(end)],'linewidth',3)
