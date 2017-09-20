function extract_bat(fn,DIAG)
bat_pos_dir='..\bat_pos\';
% DIAG=1;
%   % just testing with this file so you don't have to get a dialog box
%   fn='F:\small_space_beampattern\Test\20150820\Session 2\Trial01.c3d';

if ( nargin == 0 || isempty(fn) )
  manual_load=1;
  if ispref('small_space_beam') && ispref('small_space_beam','pname')
    pname=getpref('small_space_beam','pname');
  else
    pname='';
  end
  [fname,pname]=uigetfile('*.c3d',[],pname);
  if isequal(fname,0)
    return;
  end
  fn=[pname,fname];
  setpref('small_space_beam','pname',pname)
else
  manual_load=0;
end

[pname,fname] = fileparts(fn);
D=strsplit(pname,'\');
dateindx= ~cellfun(@isempty,strfind(D,'2015'));
datename=D{dateindx};

%grabbing band from description in enf file
enf_file=dir([pname '\' D{end} '*.enf']);
fileID = fopen([pname '\' enf_file(1).name],'r');
E=textscan(fileID,'%s','Delimiter','\n');
fclose(fileID);

species={'Rousettus','Eptesicus','Carollia','Hippo'};

rownum= ~cellfun('isempty',strfind(E{1},'DESCRIPTION'));
F=strsplit(E{1}{rownum},...
  'DESCRIPTION=');
banddata=F{2};
banddata=regexprep(banddata,species,{'','','',''},'ignorecase');
banddata=regexprep(banddata,' ','_');
if banddata(1)=='_'
  banddata(1)='';
end

sp='';
for SS=1:length(species)
  idx=strfind(F{2},species{SS});
  idx_lower=strfind(F{2},lower(species{SS}));
  if ~isempty(idx) || ~isempty(idx_lower)
    sp=[lower(species{SS}) '_'];
  end
end

if isequal(sp,'') && isequal(banddata(1:2),'LB')
  sp=['eptesicus' '_'];
end

trialnum=regexprep(fname,'Trial','');

newfname=[bat_pos_dir sp strjoin({datename,banddata,trialnum,'bat_pos'},'_') '.mat'];
newfname=regexprep(newfname,'__','_');

if exist(newfname,'file') && ~manual_load
  return
end

[point_array, frame_rate, trig_rt, trig_sig, start_f, end_f] = lc3d( fn );

if isempty(point_array)
  return;
end
point_names=cellfun(@(c) c.name,point_array,'uniformoutput',0);

%extracting
markers={'Tip','Left','Right'};
bat_pos={};
for mm=1:length(markers)
  lab=~cellfun(@isempty,strfind(point_names,markers{mm}));
  if ~isempty(find(lab, 1))
    bat_pos{mm}=point_array{lab}.traj./1e3;
    bat_pos{mm}(bat_pos{mm}==0)=nan;
  end
end

if isempty(bat_pos) || ~isempty(find(cellfun('isempty',bat_pos),1))
  return
end

%creating a plane
frames_w_alldata=find(isfinite(bat_pos{1}(:,1)) & isfinite(bat_pos{2}(:,1))...
  & isfinite(bat_pos{3}(:,1)));

flA=[0 0 0];
flB=[1 0 1];
flC=[-1 0 1];
nnfl=cross(flB-flA,flC-flA);

nn=nan(size(bat_pos{1}));
th_roll=nan(size(bat_pos{1},1),1);
head_vec=nan(size(bat_pos{1}));
for fr=frames_w_alldata'
  %three points on the triangle
  A=bat_pos{1}(fr,:);
  B=bat_pos{2}(fr,:);
  C=bat_pos{3}(fr,:);
  
  %get roll vector from plane
  nn(fr,:)=cross(B-A,C-A);
  th_roll(fr)=atan2(norm(cross(nn(fr,:),nnfl)),dot(nn(fr,:),nnfl));
  
  %direction vector
  midptBC=(B+C)./2;
  head_vec(fr,:)=A-midptBC;
end

%plotting
if nargin > 2 && DIAG && ~isempty(bat_pos)
  cols='rgb';
  eframe=frames_w_alldata(end);
  bat=bat_pos{1};
  
  figure(1), clf; set(gcf,'pos',[10 40 520 480]); hold on;
  for mm=1:length(markers)
    plot3(bat_pos{mm}(:,3),bat_pos{mm}(:,1),bat_pos{mm}(:,2),[cols(mm) '.'])
  end
  plot3(bat(:,3),bat(:,1),bat(:,2),'color',[.5 .5 .5])
  for fr=frames_w_alldata'
    plot3([bat(fr,3) bat(fr,3)+head_vec(fr,3)],...
      [bat(fr,1) bat(fr,1)+head_vec(fr,1)],...
      [bat(fr,2) bat(fr,2)+head_vec(fr,2)],'k');
    plot3([bat_pos{1}(fr,3) bat_pos{2}(fr,3)],...
      [bat_pos{1}(fr,1) bat_pos{2}(fr,1)],...
      [bat_pos{1}(fr,2) bat_pos{2}(fr,2)],'color',[.5 .5 .5])
    plot3([bat_pos{1}(fr,3) bat_pos{3}(fr,3)],...
      [bat_pos{1}(fr,1) bat_pos{3}(fr,1)],...
      [bat_pos{1}(fr,2) bat_pos{3}(fr,2)],'color',[.5 .5 .5])
  end
  
  text(bat_pos{mm}(frames_w_alldata(1),3),...
    bat_pos{mm}(frames_w_alldata(1),1),...
    bat_pos{mm}(frames_w_alldata(1),2),'Start');
  text(bat_pos{mm}(frames_w_alldata(end),3),...
    bat_pos{mm}(frames_w_alldata(end),1),...
    bat_pos{mm}(frames_w_alldata(end),2),'End');
%   plot3(flA(1),flA(2),flA(3),'or');
%   plot3(flB(1),flB(2),flB(3),'or');
%   plot3(flC(1),flC(2),flC(3),'or');

%   load_mic_pos
  load('..\mic_pos\mic_pos_20150825_S1_02&04.mat');
  plot3(mic_pos(:,3),mic_pos(:,1),mic_pos(:,2),'+k');
  axis equal, grid on, view(3)
  
  figure(2); clf; set(gcf,'pos',[10 610 520 480]);
  plot(th_roll/pi*180)
  
  figure(3); clf; set(gcf,'pos',[500 40 520 480]);
  hold on;
  plot3(A(1),A(2),A(3),'or');
  plot3(B(1),B(2),B(3),'og');
  plot3(C(1),C(2),C(3),'ob');
  plot3(midptBC(1),midptBC(2),midptBC(3),'ok');
  pl_eframe=nn(eframe,:).*20;
  plot3([A(1) A(1)+pl_eframe(1)],...
    [A(2) A(2)+pl_eframe(2)],...
    [A(3) A(3)+pl_eframe(3)],'k')
  plot3([A(1) A(1)+head_vec(fr,1)],...
    [A(2) A(2)+head_vec(fr,2)],...
    [A(3) A(3)+head_vec(fr,3)],'g');
  axis equal, grid on;
  
  plot3(flA(1),flA(2),flA(3),'or');
  plot3(flB(1),flB(2),flB(3),'or');
  plot3(flC(1),flC(2),flC(3),'or');
  axis equal, grid on;
  
  figure(4); clf; set(gcf,'pos',[500 40 520 480]);
  plot3([0 nn(eframe,1)],[0 nn(eframe,2)],[0 nn(eframe,3)],'k')
  axis equal, grid on
end

%saving
if ~isempty(bat_pos)
  if exist(newfname,'file')
    b=questdlg(['File ' newfname ' already exists. Overwrite?'], 'Overwrite?', ...
      'Yes','No','No');
    switch b
      case 'No'
        return
    end
  end
  save(newfname,'bat_pos','markers',...
    'nn','th_roll','head_vec','frame_rate')
  disp(['Saved ' newfname])
end