%merge_audio_miro
clear;

%using 15fps video

bat='36134'; %(session 2)

dates={'20150910','20150825'};
trials={'2','3'};

for dd=1:length(dates)
  for tt=1:length(trials)
    
    apath=['..\mic_recordings\' dates{dd} '_calib_mic\rousettus_36134_matfile\'];
    
    if dd==2
      insert='_';
    else
      insert='';
    end
    afname=['rousettus_' bat insert trials{tt} '_mic_data.wav'];
    
    vpath=['..\miro\' dates{dd} '\'];

    camname='13817';
    miro_file={'13479','13480';'13218','13219'};
    
    vfname=['c_' camname '_' miro_file{dd,tt} '.mp4'];
    if exist([vpath vfname],'file') && exist([apath afname],'file')
      %changes frame rate without recompressing, https://gpac.wp.mines-telecom.fr/mp4box/
      system(['F:\video_tools\MP4Box -add ' vpath vfname ...
        '#video -raw 1 -new test'])
      system(['F:\video_tools\MP4Box -add test_track1.h264:fps=15 -new '...
        vpath vfname(1:end-4) '_slow.mp4'])
      
      system(['F:\video_tools\ffmpeg_git\bin\ffmpeg -y -i ' ...
        vpath vfname(1:end-4) '_slow.mp4'...
        ' -i ' apath afname ...
        ' -c:v copy -c:a libfdk_aac -b:a 32k ' ...
        '..\merged_videos\' bat '_' dates{dd} '_' trials{tt} '_C' camname '.mp4'])
    end
    
    camname='13818';
    miro_file={'500','501';'234','235'};
    
    vfname=['c_' camname '_' miro_file{dd,tt} '.mp4'];
    if exist([vpath vfname],'file') && exist([apath afname],'file')
      %changes frame rate without recompressing, https://gpac.wp.mines-telecom.fr/mp4box/
      system(['F:\video_tools\MP4Box -add ' vpath vfname ...
        '#video -raw 1 -new test'])
      system(['F:\video_tools\MP4Box -add test_track1.h264:fps=15 -new '...
        vpath vfname(1:end-4) '_slow.mp4'])
      
      system(['F:\video_tools\ffmpeg_git\bin\ffmpeg -y -i ' ...
        vpath vfname(1:end-4) '_slow.mp4'...
        ' -i ' apath afname ...
        ' -c:v copy -c:a libfdk_aac -b:a 32k ' ...
        '..\merged_videos\' bat '_' dates{dd} '_' trials{tt} '_C' camname '.mp4'])
    end

    camname='14785';
    miro_file={'5438','5439';'5176','5177'};
    
    vfname=['c_' camname '_' miro_file{dd,tt} '.mp4'];
    if exist([vpath vfname],'file') && exist([apath afname],'file')
      %changes frame rate without recompressing, https://gpac.wp.mines-telecom.fr/mp4box/
      system(['F:\video_tools\MP4Box -add ' vpath vfname ...
        '#video -raw 1 -new test'])
      system(['F:\video_tools\MP4Box -add test_track1.h264:fps=15 -new '...
        vpath vfname(1:end-4) '_slow.mp4'])
      
      system(['F:\video_tools\ffmpeg_git\bin\ffmpeg -y -i ' ...
        vpath vfname(1:end-4) '_slow.mp4'...
        ' -i ' apath afname ...
        ' -c:v copy -c:a libfdk_aac -b:a 32k ' ...
        '..\merged_videos\' bat '_' dates{dd} '_' trials{tt} '_C' camname '.mp4'])
    end

    camname='17619';
    miro_file={'1578','1579';'1313','1314'};
    
    vfname=['c_' camname '_' miro_file{dd,tt} '.mp4'];
    if exist([vpath vfname],'file') && exist([apath afname],'file')
      %changes frame rate without recompressing, https://gpac.wp.mines-telecom.fr/mp4box/
      system(['F:\video_tools\MP4Box -add ' vpath vfname ...
        '#video -raw 1 -new test'])
      system(['F:\video_tools\MP4Box -add test_track1.h264:fps=15 -new '...
        vpath vfname(1:end-4) '_slow.mp4'])
      
      system(['F:\video_tools\ffmpeg_git\bin\ffmpeg -y -i ' ...
        vpath vfname(1:end-4) '_slow.mp4'...
        ' -i ' apath afname ...
        ' -c:v copy -c:a libfdk_aac -b:a 32k ' ...
        '..\merged_videos\' bat '_' dates{dd} '_' trials{tt} '_C' camname '.mp4'])
    end
  end
end




