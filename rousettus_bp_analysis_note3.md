# Rousettus beampattern processing and modeling ---  2017


## 2017/01/16
### Move BEM results off from server
Move most BEM calculation results off from the gal-e server. Only keep the ones used in the paper:

- **20160820_bullethead-sc-colony-0.5mm**: bullethead scaled to match dimension of Ra-colony, no rotation around y-axis
- **20160829_Ra-colony-noear-0.5mm**: Ra-colony bat head *without* ears (rotated to have mouth more parallel to x-axis)
- <span style="background-color:yellow">**20160917_Ra-colony-rotear-0.5mm**</span>: Ra-colony bat head *with ears rotated slightly downward* (rotated to have mouth more parallel to x-axis) --> this slight ear rotation eliminated the 25 kHz anomaly --> this is the final mesh used in Fig. 4 in the paper

 

## 2017/03/06
### Reprocessed Rousettus bp files
* Ben found out that the mic beampattern compensation was not done correctly: the mic_bp file has the beampattern values in degrees but in `bp_proc` the angles are in radians. Wu-Jung later found out that in `bp_proc` the function `bat2mic` was problematic that the vector `mic_vec_norm` was not normalized before used to calculate the polar angle. These are corrected as of 2017/03/01 and changed were pushed up to GitHub repo.
* Wu-Jung reprocessed Rousettus bp files on 2017/03/06. Results are stored in the following folders:
	* `proc_output_rousettus_new_raw`: raw results after `bp_proc`
	* `proc_output_rousettus_new_isgood call`: results after applying `chk_isgoodcall`
	* `proc_output_rousettus_new_checked`: results after copying the good calls from what Ben had checked before (`proc_output_rousettus_checked`). The copying was done using code `copy_isgoodcall`



## 2017/09/20
### Quantification of measured and modeled beampattern
* Beam center for individual click measurements `fig_single_bp_on_globe_20170920`:
	* if use **max point _after_ interpolation ('rx')** as the beam center: Observed that the beam centers are either on top of the mic with max receiving level or between 2-3 mics with highest receicing levels. This means the beam center found in individual measured beampatterns are strongly affected by spatial under-sampling.
	* if use **center of the best-fitting ellipse ('ro')** as the beam center, the locations of the beam center are more consistent with what one would expect by only looking at the -3 dB contours. Using this also potentially mitigates the bias from spatial under-sampling.
* Beam center for composite clicks: `fig_composite_click_avg_bp_20170920`
	* if use **max point _after_ interpolation ('rx')** as the beam center: In the plot can see the beam center moves from outer side to inner side as frequency increases from 25 kHz to 45 kHz, but from 45 kHz to 50 kHz it's a _vertical downward_ movement. This pattern is reproduced for both left and right composite clicks.
	* if use **center of the best-fitting ellipse ('ro')** as the beam center: the locations of the ellipse centers are generally fairly close to the max beam energy location
* Beam center for model results `fig_model_steer_h_bpctr`:
	* max beam energy location jumps around a bit across all steered angles, whereas center of best-fitting ellipse is very stable --> therefore use the latter in figures
* Tomorrow:
	* finish putting beam center on all figures
	* multi-freq beam center shifts for model (both individual and composite)
	* compare beam center location between data and model
	* compare el/az ratio between data and model
	* effects of missing one or more markers on the bat head



## 2017/09/21
### Cont: Quantification of measured and modeled beampattern
* Compare the max beam energy location vs averaged location for the top 1%.
	* For composite measured clicks, the averaged location of all points >-1 dB normalized beam energy ('r^') is much closer to the location of the center of the best-fitting ellipse ('ro') than just the max beam energy location ('rx').
	<img src=./img/fig_composite_click_avg_bp_20170920_batall_bin10_th0_35kHz_avg_bp_left.png width="400">	<img src=./img/fig_composite_click_avg_bp_20170920_batall_bin10_th0_35kHz_avg_bp_right.png width="400">


## TO-DO
* Use the plotting routines in `plot_indiv_click_rorate` to update those in the beampattern processing GUI, since doing conversion to map domain and plot using `contour` is much faster than calling `contourfm`



************************************************************************
## Codes used for paper figures

### Fig. 2-3
Experimental data, beampattern and multi-freq structure for both individual and composite clicks

- `fig_single_bp_on_globe_20160721`
- `fig_composite_click_avg_bp_20160808`
- `fig_el_cut_indiv_clicks_20160809` --> identical to `el_cut_indiv_clicks_20160722`, only file saving modifications
- `fig_multifreq_indiv_click_cntr_201600809`
- `fig_multifreq_composite_click_cntr_20160809`

Decided to use frequency spacing 25:5:55 kHz for multi-frequency plots. The reason is that the composite -3 dB contour for 20 kHz is very twisted for the threshold=0 and binsize=10deg case. The 20 kHz contour looks fine with higher threshold, but it's good to keep it consistent with the beampattern plot (Fig. 2A, 2B).

### Fig. 4

Code to plot mesh shapes:
- `fig_mesh_jigsaw` --> main code to plot mesh
- `draw_jigsaw_mesh` --> called by `fig_mesh_jigsaw` to draw mesh



************************************************************************
## List of BEM models calculated to date
Matching date/folder with model shape. Note all shapes are rotated such that the front faces +X axis

- **20160421**: bullethead (`bullethead_new_0.1.msh`), centered but not tilted
- **20160501**: bullethead (`bullethead_new_0.1.msh`), centered and tilted toward -Z-axis for -10 deg
- **20160727**: piston on square plate (`piston_new_0.1.msh`)
- **20160729**: circular piston shape (`piston_circ_0.05.msh`)
- **20160801**: bat head Ra\_224\_14k (`Rousettus_aegyptiacus_14k_0.88mm-remesh_fill_rot.msh`), centered and corrected for the head orientation
- **20160801b**: bat head Ra\_224\_7k (`Rousettus_aegyptiacus_7k_1.4mm-remesh_fill_rot.msh`), centered and corrected for the head orientation
- **20160803_Ra224-0.5mm**: Ra224 bat head, rotated around y-axis for -10 deg
- **20160803_Ra224-0.6mm**: Ra224 bat head, rotated around y-axis for -10 deg
- **20160811_Ra224-0.6mm**: Ra224 bat head, rotated around y-axis for -20 deg (mouth more parallel to x-axis)
- **20160812_Ra224-0.5mm**: Ra224 bat head, rotated around y-axis for -20 deg (mouth more parallel to x-axis)
- **20160812_bullethead-sc-0.5mm**: bullethead scaled to match dimension of Ra224, rotated around y-axis for 5 deg
- **20160814_bullethead-sc-0.6mm**: bullethead scaled to match dimension of Ra224, rotated around y-axis for 5 deg
- **20160814_Ra224-ear-0.5142mm**: Ra224 with ears, rotated around y-axis for -20 deg (mouth more parallel to x-axis)
- **20160817_Ra-colony-0.5mm**: Ra-colony bat head with ears (rotated to have mouth more parallel to x-axis)
- **20160817_Ra-colony-0.7mm**: Ra-colony bat head with ears (rotated to have mouth more parallel to x-axis)
- **20160817_bullethead-sc-colony-0.5mm**: bullethead scaled to match dimension of Ra-colony, rotated around y-axis for 5 deg
- **20160820_bullethead-sc-colony-0.5mm**: bullethead scaled to match dimension of Ra-colony, no rotation around y-axis
- **20160829_Ra-colony-0.7mm**: Ra-colony bat head with ears (rotated to have mouth more parallel to x-axis), this was calculated for many fine frequencies to verify if the anomaly is only at 25 kHz or also at nearby frequencies (--> the anomaly is only at 25 kHz)
- **20160829_Ra-colony-noear-0.5mm**: Ra-colony bat head *without* ears (rotated to have mouth more parallel to x-axis)
- <span style="background-color:yellow">**20160917_Ra-colony-rotear-0.5mm**</span>: Ra-colony bat head *with ears rotated slightly downward* (rotated to have mouth more parallel to x-axis) --> this slight ear rotation eliminated the 25 kHz anomaly --> this is the final mesh used in Fig. 4 in the paper
