## Code for Rousettus beampattern analysis

- `shift_rotate_bp`: shift beampattern according to the center of the best-fitting ellipse and compensate for the ellipse rotation
- `bp_fit_ellipse_azel`: called by `shift_rotate_bp` to find the best-fitting ellipse
- `rotate_all_click`: shift and rotate all clicks and save into a formatted structure (e.g., `all_click_rotated_3bats.mat`)
- `multifreq_composite_click`: load data from all rotated clicks (e.g., `all_click_rotated_3bats.mat`) and do interpolation for composite clicks
- `multifreq_composite_click_plots`: plot multi-freq composite clicks (beampattern color map and -3dB curves)
- `multifreq_contours_indiv_clicks`: plot all multi-freq -3dB contours from individual clicks overlaying on top of one another

