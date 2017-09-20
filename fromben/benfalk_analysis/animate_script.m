clear;

save_movie=1;

% bat_type='eptesicus';
bat_type='rousettus';

%extract_miro_match_beam_dirs %only needs to be run once for each trial

for freq_desired=[35 55 70]
  
  for side_view=[0 1]
    animate_beam_dirs
  end
  
  animate_beam_pattern
  merge_videos_avisynth
end