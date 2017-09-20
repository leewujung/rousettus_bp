%allows freq_desired to be passed into script as a variable
a=[];
vars=who;
clear_vars=setdiff(vars,{'freq_desired','bat_type'});
clear(clear_vars{:});

miro_view=1; %or 2

base_dir = '..\animate_beam_dirs\';

if exist('bat_type','var') && strcmp(bat_type,'rousettus')
  mic_proc_dir='..\proc_output\';
  trials=dir([mic_proc_dir 'rousettus_20150825*.mat']);
else
  mic_proc_dir='..\proc_output_eptesicus_new\';
  trials=dir([mic_proc_dir 'eptesicus_20150824_*_mic_data_bp_proc.mat']);
  checked=1;
end

if ~exist('freq_desired','var')
  freq_desired=35; %khz
end

ftypes={'c13817','c13818','c14785','c17619',...
  ['beampattern_' num2str(freq_desired)],...
  ['side_interp_' num2str(freq_desired)],['top_interp_' num2str(freq_desired)]};

for tt=1:length(trials)
  base_fn = trials(tt).name(1:end-21);
  
  if exist('checked','var') && checked
    if exist([mic_proc_dir trials(tt).name(1:end-4) '_checked.mat'],'file')
    else
      continue
    end
  end
  
  merge_files=dir([base_dir base_fn '*.mp4']);
  fnames={merge_files.name}';
  
  fname_type=nan(length(ftypes),1);
  for ft=1:length(ftypes)
    fname_type(ft)=find(~cellfun(@isempty,strfind(fnames,ftypes{ft})));
  end
  
  %generate avisynth script to combine videos
  fileID = fopen([base_dir 'loader.avs'],'w');
  
  %top view
  fprintf(fileID,...
    ['v1=crop(FFVideoSource("' fnames{fname_type(7)} '"),0,0,-60,0)\n']);
%   fprintf(fileID,...
%     'v1=Subtitle(v1,"Top",align=8,size=36,font="arial",text_color=$000000)\n');
  
  %beampattern
  fprintf(fileID,...
    ['v2=crop(FFVideoSource("' fnames{fname_type(5)}...
    '"),80,30,-60,-60)\n']);
  
  if miro_view == 1
    %miro view 1
    fprintf(fileID,...
      ['v3=crop(FFVideoSource("' fnames{fname_type(3)}...
      '").tweak(cont=1.2,bright=0,coring=true),130,300,-300,-200)\n\n']);
  else
    %miro view 2
    fprintf(fileID,...
      ['v3=crop(FFVideoSource("' fnames{fname_type(2)}...
      '").tweak(cont=2.2,bright=35,coring=true),350,70,0,0)\n\n']);
  end
    
  
  fprintf(fileID,'asp_rat = float(v3.width) / float(v3.height)\n');
  fprintf(fileID,'rat_change = float(v2.width)/ float(v3.width)\n');
  fprintf(fileID,'v8=stackvertical(v2,Spline36Resize(v3, v2.width ,  round( ( (v3.height * rat_change ) + 2 )/2  ) * 2 ) )\n\n');

  fprintf(fileID,'rat_change = float(v8.height)/ float(v1.height)\n');
  fprintf(fileID,'stackhorizontal(Spline36Resize(v1, round( ((v1.width * rat_change)+2 )/2  )*2  , v8.height ) ,v8)\n');
  
  fclose(fileID);
  
  %32 bit because of avisynth
  [status,cmdout]=system(['C:\video_tools\ffmpeg_32\bin\ffmpeg.exe -y '...
      ' -i ' base_dir 'loader.avs' ...
      ' -i ' base_dir fnames{fname_type(7)} ...
      ' -c:v libx264 -crf 20 -pix_fmt yuv420p -c:a:1 copy '...
      base_dir 'merged\' base_fn '_' num2str(freq_desired) '_miro_view_' num2str(miro_view) '.mp4']);
end


% indx_cams=nan(length(fnames),length(ftypes));
% for cc=1:length(ftypes)
%   indx_cams(:,cc) = ~cellfun('isempty',strfind(fnames,ftypes{cc}));
% end
% 
% base_file_indx = find(indx_cams(:,5));
% for ff=base_file_indx'
%   %generate avisynth script to combine videos
%   fileID = fopen([base_dir 'loader.avs'],'w');
%   fprintf(fileID,...
%     ['v1=crop(DirectShowSource("' fnames{ff} '",fps=12),0,30,-20,0)\n']);
%   fprintf(fileID,...
%     ['v2=crop(DirectShowSource("' fnames{ff+1} '",fps=12),0,40,-20,-60)\n']);
%   
%   fprintf(fileID,...
%     ['v3=crop(DirectShowSource("' fnames{ff+4} '",fps=12),120,184,0,-10)\n']);
%   fprintf(fileID,...
%     ['v4=crop(DirectShowSource("' fnames{ff+2} ...
%     '",fps=12).tweak(cont=2,bright=-8),0,0,-120,0)\n']);
%   
%   fprintf(fileID,'v5=stackvertical(v1,v2)\n');
%   fprintf(fileID,'v6=stackvertical(v3,v4)\n');
%   fprintf(fileID,'stackhorizontal(v6,v5)\n');
% 
%   fclose(fileID);
%   
%   [status,cmdout]=system(['C:\video_tools\ffmpeg_32\bin\ffmpeg.exe -y '...
%       ' -i ' base_dir 'loader.avs' ...
%       ' -i ' base_dir fnames{ff} ...
%       ' -c:v libx264 -crf 20 -c:a copy'...
%       ' ' base_dir fnames{ff}(1:end-4) '_merged.mp4']);
% end
