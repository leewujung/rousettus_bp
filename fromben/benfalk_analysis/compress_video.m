%function status=compress_video(fname,deleteOrig,mv_file,mv_dir)
function status=compress_video(fname,deleteOrig,mv_file,mv_dir)

% system(['C:\video_tools\ffmpeg_64\bin\ffmpeg.exe -y -i "' fname ...
%   '" -c:v libxvid -q:v 4 "' fname(1:end-3) 'mkv"']);
[status,~]=system(['C:\video_tools\ffmpeg_64\bin\ffmpeg.exe -y -i "' fname ...
  '" -pix_fmt yuv420p -c:v libx264 -crf 20 "' fname(1:end-3) 'mkv"']);
if status~=0
  disp('something went wrong during compression')
  return;
end
  
if status==0 && nargin > 1
  if nargin>2 && mv_file
    movefile([fname(1:end-3) 'mkv'],mv_dir);
  end
end

if deleteOrig
  delete(fname)
end